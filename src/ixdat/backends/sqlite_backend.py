"""This module implements an SQLite database backend for ixdat

Saves and loads ixdat objects using Python's built-in ``sqlite3`` module, so no
extra dependency is needed. The table layout isn't written by hand - it's built
from the same table info each Saveable class already carries (see
:module:`~ixdat.backends.relational`). Tables are created the first time they're
needed, so any Saveable class, including ones from external plugins, works here
without extra setup.

Basic usage::

    from ixdat.db import change_database

    change_database("sqlite", db_path="my_project.sqlite")
    measurement.save()
    ...
    loaded = Measurement.get(1)
    loaded = Measurement.load("my measurement")  # or, by name

Numeric arrays (numpy data) are stored as blobs and are only read from the
database when ``.data`` is actually accessed, so opening an object stays fast
even in a large project.

Saving an object and everything it references happens in one transaction, so a
crash partway through can't leave the database half-written. The first time a
table is used, it's checked against what the current Saveable class expects;
missing columns are added automatically, and anything that can't be changed
safely raises a :class:`~ixdat.exceptions.DataBaseError` instead of silently
doing the wrong thing.

Note: one SQLiteBackend belongs to the thread that created it - this is just how
Python's sqlite3 module works. Make a separate backend instance per thread.
"""

import json
import sqlite3
from io import BytesIO
from pathlib import Path

import numpy as np

from . import relational
from .backend_base import BackendBase
from ..config import config, prompt_for_permission
from ..exceptions import DataBaseError


DATABASE_FILE_SUFFIX = ".sqlite"
SCHEMA_VERSION = 1
METADATA_TABLE = "_ixdat_metadata"


def _quote_identifier(identifier):
    """Return an SQLite identifier with embedded quotes escaped."""
    return '"' + identifier.replace('"', '""') + '"'


class SQLiteBackend(BackendBase):
    """A database backend which saves and loads ixdat objects in an SQLite file"""

    backend_type = "sqlite"

    def __init__(self, db_path=None, directory=None, project_name=None):
        """Initiate the backend, connecting to (and if needed creating) the database

        Args:
            db_path (Path or str): The path to the SQLite database file. By
                default it is named after the project and put in ixdat's
                standard data directory, mirroring the directory backend.
            directory (Path): The directory for the default database file
            project_name (str): The project name for the default database file
        """
        self._is_memory = str(db_path) == ":memory:"
        if self._is_memory:
            self.db_path = ":memory:"
        elif db_path:
            self.db_path = Path(db_path)
        else:
            directory = Path(directory or config.standard_data_directory)
            project_name = project_name or config.default_project_name
            self.db_path = directory / (project_name + DATABASE_FILE_SUFFIX)
        self._resolved_path = None  # what identifies this database, see __eq__
        if not self._is_memory:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            # resolved once here: two backends on one file must compare and hash
            # alike however their paths were written, and resolving hits the disk
            self._resolved_path = self.db_path.resolve()
        self.connection = sqlite3.connect(self.db_path)
        self.connection.execute("PRAGMA foreign_keys = ON")
        self._closed = False
        self._ensured_classes = set()
        self._ensured_schema_definitions = set()
        try:
            self._initialize_schema_metadata()
        except Exception:
            self.close()
            raise
        super().__init__()

    @property
    def address(self):
        """The path to the SQLite database file"""
        return str(self.db_path)

    def close(self):
        """Close the connection to the database file"""
        if not self._closed:
            self.connection.close()
            self._closed = True

    def __enter__(self):
        """Return this backend for use as a context manager."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close the connection when leaving a context manager."""
        self.close()

    def __eq__(self, other):
        """Two SQLiteBackends are equivalent if they use the same database file"""
        if other is self:
            return True
        if other.__class__ is not self.__class__:
            return False
        if self._is_memory or other._is_memory:
            # Separate ":memory:" connections are separate databases.
            return False
        return other._resolved_path == self._resolved_path

    def __hash__(self):
        """Hash by the database file, so that equal backends hash alike

        Defining __eq__ sets __hash__ to None unless it is defined too, which would
        make this backend unhashable. Saveable does the same thing for the same
        reason. Two ":memory:" backends are never equal, so they hash by identity.
        """
        if self._is_memory:
            return object.__hash__(self)
        return hash(self._resolved_path)

    # ------- saving  ------- #

    def save(self, obj, force=False, no_updates=True):
        """Atomically save an object and all objects it references.

        Args:
            obj (Saveable): the object to save
            force (bool): Whether to force updates if the object is already saved
            no_updates (bool): Whether to allow updates if the object is already
                saved. If both force and no_updates are False, the user will be
                prompted on whether to save.
        """
        plan = self._build_save_plan(obj, force=force, no_updates=no_updates)
        for planned_obj, action in plan:
            if action != "skip":
                self._ensure_tables(type(planned_obj))

        original_states = []
        root_result = None
        try:
            # Children come before their parents in "plan", so by the time we
            # insert a row, anything it points to already has its id. Doing all
            # inserts/updates on one connection makes them one transaction: if
            # anything below fails, SQLite rolls every row in "plan" back.
            with self.connection as connection:
                for planned_obj, action in plan:
                    if action == "skip":
                        continue
                    if action == "update":
                        self._update_row(connection, planned_obj)
                        result = planned_obj.id
                    else:
                        original_states.append(
                            (planned_obj, planned_obj._backend, planned_obj._id)
                        )
                        result = self._add_row(connection, planned_obj)
                        planned_obj.set_id(result)
                        planned_obj.set_backend(self)
                    if planned_obj is obj:
                        root_result = result
        except Exception:
            # The database rollback undoes the inserted rows; undo the matching
            # id/backend changes on the in-memory objects too, so they don't end
            # up claiming to be saved when they aren't.
            for saved_obj, old_backend, old_id in reversed(original_states):
                saved_obj._backend = old_backend
                saved_obj._id = old_id
            raise
        return root_result

    def _build_save_plan(self, root, force, no_updates):
        """Return unique objects in child-first order with insert/update actions."""
        plan = []
        visited = set()  # objects already added to "plan"
        visiting = set()  # objects currently being visited, i.e. above us in the walk
        existing_tables = self._existing_tables()

        def visit(obj, updates_forbidden):
            object_key = id(obj)
            if object_key in visited:
                return  # already planned (e.g. two parents sharing one child)
            if object_key in visiting:
                # we're already in the middle of visiting this object further up
                # the call stack, i.e. it refers back to itself somewhere
                raise DataBaseError(
                    f"Can't save a cyclic object graph rooted at {root!r}."
                )
            visiting.add(object_key)
            # visit every child first, so they end up earlier in "plan"
            if obj.child_attrs:
                for child_list_name in obj.child_attrs:
                    for child_obj in getattr(obj, child_list_name) or []:
                        visit(child_obj, updates_forbidden=True)
            visiting.remove(object_key)
            visited.add(object_key)
            plan.append(
                (
                    obj,
                    self._save_action(obj, force, updates_forbidden, existing_tables),
                )
            )

        visit(root, updates_forbidden=no_updates)
        return plan

    def _save_action(self, obj, force, updates_forbidden, existing_tables):
        """Return ``insert``, ``update``, or ``skip`` for one planned object."""
        row_exists = False
        if obj.backend is self and obj.table_name in existing_tables:
            row_exists = self.connection.execute(
                f"SELECT 1 FROM {_quote_identifier(obj.table_name)} " 'WHERE "id" = ?',
                (obj.id,),
            ).fetchone()
        if row_exists:
            if force:
                return "update"
            if not updates_forbidden and prompt_for_permission(
                f"Are you sure you would like to overwrite "
                f"{self} table={obj.table_name} id={obj.id} with {obj}? "
                f"(You can use save() with force=True to suppress this.)"
            ):
                return "update"
            return "skip"
        return "insert"

    def _add_row(self, connection, obj):
        """Insert one object's rows using an existing transaction."""
        main_schema = relational.main_table_schema(type(obj))
        main_dict = obj.get_main_dict()
        columns = main_schema.data_columns
        values = [self._encode(column, main_dict[column.name]) for column in columns]
        cursor = connection.execute(
            _insert_sql(main_schema.name, [column.name for column in columns]),
            values,
        )
        i = cursor.lastrowid
        self._insert_extras(connection, obj, i)
        return i

    def _update_row(self, connection, obj):
        """Update one object's rows using an existing transaction."""
        main_schema = relational.main_table_schema(type(obj))
        main_dict = obj.get_main_dict()
        columns = main_schema.data_columns
        values = [self._encode(column, main_dict[column.name]) for column in columns]
        assignments = ", ".join(
            f"{_quote_identifier(column.name)} = ?" for column in columns
        )
        connection.execute(
            f"UPDATE {_quote_identifier(main_schema.name)} "
            f'SET {assignments} WHERE "id" = ?',
            values + [obj.id],
        )
        for schema in relational.extension_table_schemas(type(obj)):
            connection.execute(
                f'DELETE FROM {_quote_identifier(schema.name)} WHERE "id" = ?',
                (obj.id,),
            )
        for linker in relational.linker_table_schemas(type(obj)):
            connection.execute(
                f"DELETE FROM {_quote_identifier(linker.name)} WHERE "
                f"{_quote_identifier(linker.owner_column)} = ?",
                (obj.id,),
            )
        self._insert_extras(connection, obj, obj.id)

    def _insert_extras(self, connection, obj, i):
        """Insert the object's rows in its extension and linker tables"""
        # extension tables hold the extra columns of a subclass, one row per object,
        # keyed by the same id as the main table row:
        for schema in relational.extension_table_schemas(type(obj)):
            columns = schema.data_columns
            values = [
                self._encode(column, getattr(obj, column.name)) for column in columns
            ]
            connection.execute(
                _insert_sql(schema.name, ["id"] + [column.name for column in columns]),
                [i] + values,
            )
        # linker tables hold one row per reference to another object, e.g. one row
        # per data series a measurement owns, in the order they should be read back:
        for linker in relational.linker_table_schemas(type(obj)):
            identities = getattr(obj, linker.id_attr)
            if identities is None:
                continue
            if not linker.is_list:
                identities = [identities]
            for position, identity in enumerate(identities):
                linked_id = self._dereference(identity)
                if not isinstance(linked_id, int):
                    raise DataBaseError(
                        f"{obj!r} refers to {linker.id_attr}={identity}, which is "
                        f"not the id of a row saved in {self}. Do the referenced "
                        "objects live in another backend?"
                    )
                connection.execute(
                    _insert_sql(
                        linker.name,
                        [linker.owner_column, "position", linker.linked_column],
                    ),
                    (i, position, linked_id),
                )

    # ------- loading  ------- #

    def get(self, cls, i):
        """Return the object of Saveable class cls built from the rows with id=i

        We don't know the row's exact subclass before reading it (an id=5 in
        "measurement" could be a plain ECMeasurement or an ECMSMeasurement), so we
        check the extension/linker tables of every class that shares cls's main
        table, not just cls's own. Only the tables the row's real class wrote to
        will actually have a row for this id; the rest are skipped. Numerical data
        isn't loaded here - see `load_obj_data`.
        """
        self._ensure_tables(cls)
        main_schema = relational.main_table_schema(cls)
        extension_schemas, linker_schemas = relational.family_table_schemas(cls)
        existing_tables = self._ensure_existing_family_tables(
            extension_schemas + linker_schemas
        )
        # Leave NDARRAY columns out of the query; `load_obj_data` gets them lazily:
        columns = [
            column for column in main_schema.data_columns if column.dtype != "NDARRAY"
        ]
        row = self.connection.execute(
            _select_sql(main_schema.name, [column.name for column in columns])
            + ' WHERE "id" = ?',
            (i,),
        ).fetchone()
        if row is None:
            raise DataBaseError(
                f"{self} has no row with id={i} in table '{main_schema.name}'"
            )
        obj_as_dict = {
            column.name: self._decode(column, value)
            for column, value in zip(columns, row)
        }
        for column in main_schema.data_columns:
            if column.dtype == "NDARRAY":
                obj_as_dict[column.name] = None  # signals lazy loading
        for schema in extension_schemas:
            if schema.name not in existing_tables:
                continue
            columns = schema.data_columns
            row = self.connection.execute(
                _select_sql(schema.name, [column.name for column in columns])
                + ' WHERE "id" = ?',
                (i,),
            ).fetchone()
            if row is not None:
                obj_as_dict.update(
                    {
                        column.name: self._decode(column, value)
                        for column, value in zip(columns, row)
                    }
                )
        for linker in linker_schemas:
            if linker.name not in existing_tables:
                continue
            rows = self.connection.execute(
                _select_sql(linker.name, [linker.linked_column])
                + f" WHERE {_quote_identifier(linker.owner_column)} = ? "
                'ORDER BY "position"',
                (i,),
            ).fetchall()
            if rows:
                linked_ids = [row[0] for row in rows]
                obj_as_dict[linker.id_attr] = (
                    linked_ids if linker.is_list else linked_ids[0]
                )
        obj = cls.from_dict(obj_as_dict)
        obj.set_backend(self)
        obj.set_id(i)
        return obj

    def load(self, cls, name):
        """Return the most recently saved object of Saveable class cls with the name"""
        self._ensure_tables(cls)
        row = self.connection.execute(
            f'SELECT "id" FROM {_quote_identifier(cls.table_name)} '
            'WHERE "name" = ? '
            'ORDER BY "id" DESC LIMIT 1',
            (name,),
        ).fetchone()
        if row is None:
            raise DataBaseError(
                f"{self} has no row named '{name}' in table '{cls.table_name}'"
            )
        return self.get(cls, row[0])

    def load_obj_data(self, obj):
        """Return the numerical data of an object, from its NDARRAY column"""
        main_schema = relational.main_table_schema(type(obj))
        for column in main_schema.data_columns:
            if column.dtype == "NDARRAY":
                row = self.connection.execute(
                    _select_sql(main_schema.name, [column.name]) + ' WHERE "id" = ?',
                    (obj.id,),
                ).fetchone()
                return self._decode(column, row[0]) if row else None

    def contains(self, table_name, i):
        """Check if id `i` is already a principle key in the table named `table_name`"""
        if table_name not in self._existing_tables():
            return False
        row = self.connection.execute(
            f'SELECT 1 FROM {_quote_identifier(table_name)} WHERE "id" = ?', (i,)
        ).fetchone()
        return row is not None

    # ------- schema  ------- #

    def _initialize_schema_metadata(self):
        """Create and validate the small, backend-owned schema metadata table."""
        metadata_table = _quote_identifier(METADATA_TABLE)
        with self.connection as connection:
            connection.execute(
                f"CREATE TABLE IF NOT EXISTS {metadata_table} ("
                '"key" TEXT PRIMARY KEY, "value" TEXT NOT NULL)'
            )
            row = connection.execute(
                f'SELECT "value" FROM {metadata_table} ' 'WHERE "key" = ?',
                ("schema_version",),
            ).fetchone()
            if row is None:
                # no version stored yet: an older database, or a brand new one.
                # Either way, _ensure_tables() will check/fix each table as it's used.
                connection.execute(
                    f'INSERT INTO {metadata_table} ("key", "value") ' "VALUES (?, ?)",
                    ("schema_version", str(SCHEMA_VERSION)),
                )
                return
            try:
                database_version = int(row[0])
            except (TypeError, ValueError):
                raise DataBaseError(
                    f"Invalid ixdat schema version {row[0]!r} in {self.db_path}."
                )
            if database_version != SCHEMA_VERSION:
                relation = "newer" if database_version > SCHEMA_VERSION else "older"
                raise DataBaseError(
                    f"The database {self.db_path} uses {relation} ixdat schema "
                    f"version {database_version}; this ixdat supports version "
                    f"{SCHEMA_VERSION}. Migrate the database before opening it."
                )

    def _ensure_tables(self, cls):
        """Create, validate, and additively migrate tables needed by ``cls``."""
        # this only needs doing once per class per backend instance, since a table
        # doesn't change shape again while this backend stays open:
        if cls in self._ensured_classes:
            return
        existing_tables = self._existing_tables()
        ensured_definitions = []
        with self.connection as connection:
            for schema in relational.table_schemas_of(cls):
                if schema.name in existing_tables:
                    self._validate_or_migrate_table(connection, schema)
                else:
                    connection.execute(_create_table_sql(schema))
                    existing_tables.add(schema.name)
                self._ensure_indexes(connection, schema)
                ensured_definitions.append((schema.name, _create_table_sql(schema)))
        self._ensured_schema_definitions.update(ensured_definitions)
        self._ensured_classes.add(cls)

    def _ensure_existing_family_tables(self, schemas):
        """Validate subclass tables and return the current table-name set.

        Called from `get()`, which - unlike `_ensure_tables()` - needs to check
        tables belonging to *other* classes than the one it was asked for (see the
        comment in `get()`). We track each table by its exact expected DDL here
        rather than by class, since several classes can share one such table.
        """
        existing_tables = self._existing_tables()
        ensured_definitions = []
        with self.connection as connection:
            for schema in schemas:
                schema_definition = (schema.name, _create_table_sql(schema))
                if (
                    schema.name not in existing_tables
                    or schema_definition in self._ensured_schema_definitions
                ):
                    continue  # table doesn't exist yet, or was already checked
                self._validate_or_migrate_table(connection, schema)
                self._ensure_indexes(connection, schema)
                ensured_definitions.append(schema_definition)
        self._ensured_schema_definitions.update(ensured_definitions)
        return existing_tables

    def _validate_or_migrate_table(self, connection, schema):
        """Add missing nullable columns and reject incompatible table layouts."""
        table_info = {
            row[1]: {"type": row[2].upper(), "not_null": row[3], "pk": row[5]}
            for row in connection.execute(
                f"PRAGMA table_info({_quote_identifier(schema.name)})"
            ).fetchall()
        }
        if isinstance(schema, relational.LinkerTableSchema):
            # linker tables aren't migrated column-by-column like regular tables
            # below: if their shape is wrong at all, we just ask for a manual fix.
            expected_names = {
                schema.owner_column,
                "position",
                schema.linked_column,
            }
            if not expected_names.issubset(table_info):
                self._incompatible_schema(
                    schema.name,
                    "linker columns are missing or use obsolete names",
                )
            for column_name in expected_names:
                column = table_info[column_name]
                if column["type"] != "INTEGER" or not column["not_null"]:
                    self._incompatible_schema(
                        schema.name,
                        f"linker column {column_name!r} must be a non-null INTEGER",
                    )
            expected_pk = [schema.owner_column, "position"]
            actual_pk = [
                name
                for name, info in sorted(
                    table_info.items(), key=lambda item: item[1]["pk"] or 999
                )
                if info["pk"]
            ]
            if actual_pk != expected_pk:
                self._incompatible_schema(
                    schema.name,
                    f"primary key is {actual_pk}, expected {expected_pk}",
                )
            self._validate_foreign_key(
                connection,
                schema.name,
                schema.owner_column,
                schema.owner_table,
                "id",
            )
            self._validate_foreign_key(
                connection,
                schema.name,
                schema.linked_column,
                schema.linked_table,
                "id",
            )
            return

        # regular (main/extension) table: check each expected column one by one.
        for column in schema.columns:
            if column.name not in table_info:
                if column.name == "id":
                    self._incompatible_schema(
                        schema.name, 'primary-key column "id" missing'
                    )
                # a new column is always safe to add: it's nullable, so existing
                # rows just get NULL for it - nothing to fill in or guess at.
                connection.execute(
                    f"ALTER TABLE {_quote_identifier(schema.name)} ADD COLUMN "
                    + _column_definition(column, include_primary_key=False)
                )
                continue
            actual = table_info[column.name]
            expected_type = column.dtype or ""
            if expected_type and actual["type"] != expected_type:
                self._incompatible_schema(
                    schema.name,
                    f"column {column.name!r} has type {actual['type']!r}, "
                    f"expected {expected_type!r}",
                )
            if column.name == "id" and not actual["pk"]:
                self._incompatible_schema(schema.name, 'column "id" is not primary key')
            if column.foreign_key:
                self._validate_foreign_key(
                    connection, schema.name, column.name, *column.foreign_key
                )

    def _validate_foreign_key(
        self, connection, table_name, column_name, foreign_table, foreign_column
    ):
        """Require one expected foreign-key relationship on an existing table."""
        foreign_keys = connection.execute(
            f"PRAGMA foreign_key_list({_quote_identifier(table_name)})"
        ).fetchall()
        relationship = (column_name, foreign_table, foreign_column)
        actual_relationships = {(row[3], row[2], row[4]) for row in foreign_keys}
        if relationship not in actual_relationships:
            self._incompatible_schema(
                table_name,
                f"foreign key {column_name!r} -> "
                f"{foreign_table!r}.{foreign_column!r} is missing",
            )

    def _incompatible_schema(self, table_name, reason):
        """Raise an actionable error for a migration that cannot be done safely."""
        raise DataBaseError(
            f"Incompatible schema for table {table_name!r} in {self.db_path}: "
            f"{reason}. A manual database migration is required."
        )

    @staticmethod
    def _ensure_indexes(connection, schema):
        """Create indexes used by name lookup, dispatch, and reverse relationships."""
        if isinstance(schema, relational.LinkerTableSchema):
            indexed_columns = [(schema.linked_column,)]
        else:
            column_names = set(schema.column_names)
            indexed_columns = []
            if "name" in column_names:
                indexed_columns.append(("name", "id"))
            for discriminator in ("technique", "series_type", "calculator_type"):
                if discriminator in column_names:
                    indexed_columns.append((discriminator,))
        for columns in indexed_columns:
            index_name = "ixdat_" + schema.name + "_" + "_".join(columns)
            column_sql = ", ".join(_quote_identifier(column) for column in columns)
            connection.execute(
                f"CREATE INDEX IF NOT EXISTS {_quote_identifier(index_name)} "
                f"ON {_quote_identifier(schema.name)} ({column_sql})"
            )

    def _existing_tables(self):
        """Return the set of names of the tables in the database"""
        rows = self.connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        ).fetchall()
        return {row[0] for row in rows}

    @staticmethod
    def schema_report():
        """Return the CREATE TABLE statements of all imported Saveable classes

        This is the full relational schema of ixdat's data model, derived from
        the classes' table metadata - whether or not the tables have (yet) been
        created in this database.

        A class whose own table metadata is broken is left out of the report,
        which stays readable for every other class. Its error is raised when that
        class is itself saved or loaded.
        """
        statements = {}
        for cls in relational.saveable_classes():
            if not cls.table_name:
                continue
            try:
                schemas = relational.table_schemas_of(cls)
            except DataBaseError:
                continue
            for schema in schemas:
                statements.setdefault(schema.name, _create_table_sql(schema))
        return ";\n\n".join(statements.values()) + ";"

    # ------- value encoding  ------- #

    def _encode(self, column, value):
        """Return the database representation of a value, per its column's dtype"""
        value = self._dereference(value)
        if value is None:
            return None
        if column.dtype == "NDARRAY":
            buffer = BytesIO()
            np.save(buffer, np.asarray(value), allow_pickle=False)
            return sqlite3.Binary(buffer.getvalue())
        if column.dtype == "JSON":
            # metadata dicts sometimes contain numpy numbers/arrays, which plain
            # json.dumps can't handle - convert those to normal Python values first:
            from ..tools import to_jsonable

            return json.dumps(to_jsonable(value))
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (dict, list, tuple, np.ndarray)):
            raise DataBaseError(
                f"Can't save value {value!r} in the dynamically typed column "
                f"'{column.name}'. If this column should hold dicts/lists or numpy "
                f"arrays, give it the type 'JSON' or 'NDARRAY', respectively, in "
                "the `column_types` of the class which defines the column."
            )
        return value

    def _decode(self, column, value):
        """Return the python representation of a value, per its column's dtype"""
        if value is None:
            return None
        if column.dtype == "NDARRAY":
            return np.load(BytesIO(value), allow_pickle=False)
        if column.dtype == "JSON":
            return json.loads(value)
        return value

    def _dereference(self, value):
        """Turn a reference to another object into a plain id, for storing as a FK

        An id-holding attribute (e.g. `s_ids`) can hold either a plain int id, or,
        if the referenced object isn't in this backend yet, a `(backend, id)` pair
        (`Saveable.short_identity`). Only the plain-int case can actually be saved
        as a foreign key here, so this unwraps that pair and checks it points here.
        """
        if (
            isinstance(value, tuple)
            and len(value) == 2
            and isinstance(value[0], BackendBase)
        ):
            backend, i = value
            if backend is self or backend == self:
                return i
            raise DataBaseError(
                f"Can't save a reference to id={i} of {backend} in {self}. "
                "Save the referenced object here first."
            )
        return value


def _create_table_sql(schema):
    """Return the CREATE TABLE statement for a table or linker table schema"""
    if isinstance(schema, relational.LinkerTableSchema):
        return (
            f"CREATE TABLE IF NOT EXISTS {_quote_identifier(schema.name)} (\n"
            f"    {_quote_identifier(schema.owner_column)} INTEGER NOT NULL "
            f'REFERENCES {_quote_identifier(schema.owner_table)}("id"),\n'
            f'    "position" INTEGER NOT NULL,\n'
            f"    {_quote_identifier(schema.linked_column)} INTEGER NOT NULL "
            f'REFERENCES {_quote_identifier(schema.linked_table)}("id"),\n'
            f"    PRIMARY KEY ({_quote_identifier(schema.owner_column)}, "
            '"position")\n'
            ")"
        )
    column_definitions = [
        "    " + _column_definition(column) for column in schema.columns
    ]
    return (
        f"CREATE TABLE IF NOT EXISTS {_quote_identifier(schema.name)} (\n"
        + ",\n".join(column_definitions)
        + "\n)"
    )


def _column_definition(column, include_primary_key=True):
    """Return the DDL fragment defining one regular table column."""
    definition = _quote_identifier(column.name)
    if column.dtype:
        definition += " " + column.dtype
    if include_primary_key and column.name == "id":
        definition += " PRIMARY KEY"
    if column.foreign_key:
        foreign_table, foreign_column = column.foreign_key
        definition += (
            f" REFERENCES {_quote_identifier(foreign_table)}"
            f"({_quote_identifier(foreign_column)})"
        )
    return definition


def _insert_sql(table_name, column_names):
    """Return a parametrized INSERT statement for the named table and columns"""
    columns = ", ".join(_quote_identifier(name) for name in column_names)
    placeholders = ", ".join("?" for _ in column_names)
    return (
        f"INSERT INTO {_quote_identifier(table_name)} ({columns}) "
        f"VALUES ({placeholders})"
    )


def _select_sql(table_name, column_names):
    """Return a SELECT statement for the named table and columns"""
    columns = ", ".join(_quote_identifier(name) for name in column_names)
    return f"SELECT {columns} FROM {_quote_identifier(table_name)}"
