"""SQLite storage backend for ixdat.

SQLiteBackend stores complete ixdat object graphs in one SQLite file using the
table descriptions from ixdat.backends.relational.

See :ref:`backend` for usage, persistence guarantees, schema updates, and
connection lifetimes.
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
METADATA_TABLE_NAME = "_ixdat_metadata"
SQLITE_TYPES = {
    int: "INTEGER",
    float: "REAL",
    str: "TEXT",
    dict: "JSON",
    list: "JSON",
    tuple: "JSON",
    np.ndarray: "NDARRAY",
}


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
        # SQLite gives each ":memory:" connection its own temporary database. Give
        # each one a unique ixdat address so their saved-object identities stay apart.
        self._memory_address = f":memory:{id(self)}" if self._is_memory else None
        if self._is_memory:
            self.db_path = ":memory:"
            database_is_new = True
        elif db_path:
            self.db_path = Path(db_path)
            database_is_new = (
                not self.db_path.exists() or self.db_path.stat().st_size == 0
            )
        else:
            directory = Path(directory or config.standard_data_directory)
            project_name = project_name or config.default_project_name
            self.db_path = directory / (project_name + DATABASE_FILE_SUFFIX)
            database_is_new = (
                not self.db_path.exists() or self.db_path.stat().st_size == 0
            )
        # The resolved path identifies file storage in shares_storage_with().
        self._resolved_path = None
        if not self._is_memory:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            # Resolve once so differently written paths to one file identify the
            # same storage without repeatedly accessing the filesystem.
            self._resolved_path = self.db_path.resolve()
        self.connection = sqlite3.connect(self.db_path)
        # SQLite configures foreign-key checking per connection and leaves it off
        # by default. Enable it here so every saved id must point to an existing row:
        # https://www.sqlite.org/pragma.html#pragma_foreign_keys
        self.connection.execute("PRAGMA foreign_keys = ON")
        self._closed = False
        self._ensured_classes = set()
        self._ensured_schema_definitions = set()
        try:
            self._initialize_schema_metadata(database_is_new)
        except Exception:
            self.close()
            raise
        super().__init__()

    @property
    def address(self):
        """The file path, or the unique identity of an in-memory database."""
        return self._memory_address if self._is_memory else str(self._resolved_path)

    def close(self):
        """Close the connection to the database file"""
        if not self._closed:
            self.connection.close()
            self._closed = True

    def shares_storage_with(self, other):
        """Return whether ``other`` reaches the same SQLite database."""
        if other is self:
            return True
        if other.__class__ is not self.__class__:
            return False
        if self._is_memory or other._is_memory:
            # Separate ":memory:" connections are separate databases.
            return False
        return other._resolved_path == self._resolved_path

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
                self._ensure_tables_for_save(type(planned_obj))

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
            # Python's built-in id() identifies this exact in-memory object. That
            # lets the walk detect a cycle even when two ixdat objects have equal
            # saved values or database ids.
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
            # Visit every related object first, so it appears earlier in the plan and
            # has an id by the time this object stores the relationship.
            for related_obj in obj.iter_related_objects():
                visit(related_obj, updates_forbidden=True)
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
        # Storage identity tells us where to look. This query separately confirms
        # that the object's exact row is present through this connection.
        if self.shares_storage_with(obj.backend) and obj.table_name in existing_tables:
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
        existing_tables = self._existing_tables()
        for linker in relational.linker_table_schemas(type(obj)):
            # A connection table can be created before its linked table is needed.
            # Until that linked table exists, this object cannot have any rows in the
            # connection table to replace.
            if linker.linked_table not in existing_tables:
                continue
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
            if not linker.many:
                identities = [identities]
            for position, identity in enumerate(identities):
                linked_id = self._dereference(identity)
                if not isinstance(linked_id, int):
                    raise DataBaseError(
                        f"{obj!r} refers to {linker.id_attr}={identity}, which is "
                        f"not the id of a row saved in {self}. Do the referenced "
                        "objects live in another backend?"
                    )
                # `connection` is the sqlite3 connection opened by this backend.
                # For a measurement linked to a data series, `_insert_sql()` makes:
                # INSERT INTO "measurement_series"
                #     ("measurement_id", "position", "data_series_id")
                #     VALUES (?, ?, ?)
                # The values below give SQLite the measurement id, the series's
                # place in the measurement's list, and the data-series id.
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
        will actually have a row for this id; the rest are skipped. The ``data``
        column is loaded when the object's ``.data`` is first used - see
        `load_obj_data`.
        """
        main_schema = relational.main_table_schema(cls)
        extension_schemas, linker_schemas = relational.family_table_schemas(cls)
        existing_tables = self._validate_existing_tables(
            [main_schema] + extension_schemas + linker_schemas
        )
        if main_schema.name not in existing_tables:
            raise DataBaseError(
                f"{self} has no table for objects of type {cls.__name__}."
            )
        # `data` has ixdat's lazy-loading property. Other NumPy-array columns load
        # here with the rest of the object's values.
        columns = [
            column for column in main_schema.data_columns if column.name != "data"
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
            if column.name == "data":
                obj_as_dict["data"] = None  # signals lazy loading
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
                    linked_ids if linker.many else linked_ids[0]
                )
        obj = cls.from_dict(obj_as_dict)
        obj.set_backend(self)
        obj.set_id(i)
        return obj

    def load(self, cls, name):
        """Return the most recently saved object of Saveable class cls with the name"""
        main_schema = relational.main_table_schema(cls)
        existing_tables = self._validate_existing_tables([main_schema])
        if main_schema.name not in existing_tables:
            raise DataBaseError(
                f"{self} has no table for objects of type {cls.__name__}."
            )
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
        """Return the value stored in an object's lazy ``data`` column."""
        main_schema = relational.main_table_schema(type(obj))
        data_column = next(
            (column for column in main_schema.data_columns if column.name == "data"),
            None,
        )
        if data_column is None:
            return None
        row = self.connection.execute(
            _select_sql(main_schema.name, ["data"]) + ' WHERE "id" = ?',
            (obj.id,),
        ).fetchone()
        return self._decode(data_column, row[0]) if row else None

    def contains(self, table_name, i):
        """Check if id `i` is already a principle key in the table named `table_name`"""
        if table_name not in self._existing_tables():
            return False
        row = self.connection.execute(
            f'SELECT 1 FROM {_quote_identifier(table_name)} WHERE "id" = ?', (i,)
        ).fetchone()
        return row is not None

    # ------- schema  ------- #

    def _initialize_schema_metadata(self, database_is_new):
        """Record a new file's schema version or validate an existing file's."""
        metadata_table = _quote_identifier(METADATA_TABLE_NAME)
        if database_is_new:
            with self.connection as connection:
                connection.execute(
                    f"CREATE TABLE {metadata_table} ("
                    '"key" TEXT PRIMARY KEY, "value" TEXT NOT NULL)'
                )
                connection.execute(
                    f'INSERT INTO {metadata_table} ("key", "value") VALUES (?, ?)',
                    ("schema_version", str(SCHEMA_VERSION)),
                )
            return

        if METADATA_TABLE_NAME not in self._existing_tables():
            raise DataBaseError(
                f"The existing database {self.db_path} has no ixdat schema version. "
                "Copy its objects into a new database with an explicit migration."
            )
        try:
            row = self.connection.execute(
                f'SELECT "value" FROM {metadata_table} ' 'WHERE "key" = ?',
                ("schema_version",),
            ).fetchone()
        except sqlite3.DatabaseError as error:
            raise DataBaseError(
                f"The database {self.db_path} has invalid ixdat schema metadata."
            ) from error
        if row is None:
            raise DataBaseError(
                f"The database {self.db_path} has no ixdat schema version. "
                "Copy its objects into a new database with an explicit migration."
            )
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
                f"{SCHEMA_VERSION}. Copy its objects into a new database using "
                "a compatible ixdat version."
            )

    def _ensure_tables_for_save(self, cls):
        """Create or validate the current-version tables needed to save ``cls``."""
        # this only needs doing once per class per backend instance, since a table
        # doesn't change shape again while this backend stays open:
        if cls in self._ensured_classes:
            return
        existing_tables = self._existing_tables()
        ensured_definitions = []
        with self.connection as connection:
            for schema in relational.table_schemas_of(cls):
                if schema.name in existing_tables:
                    self._validate_table(connection, schema)
                else:
                    connection.execute(_create_table_sql(schema))
                    existing_tables.add(schema.name)
                self._ensure_indexes(connection, schema)
                ensured_definitions.append((schema.name, _create_table_sql(schema)))
        self._ensured_schema_definitions.update(ensured_definitions)
        self._ensured_classes.add(cls)

    def _validate_existing_tables(self, schemas):
        """Validate existing tables without changing the database.

        ``get()`` checks tables belonging to sibling classes before it knows the
        exact class of a saved row. Missing sibling tables are fine. The caller
        decides whether a missing main table means that the requested object type
        has never been saved.
        """
        existing_tables = self._existing_tables()
        ensured_definitions = []
        for schema in schemas:
            schema_definition = (schema.name, _create_table_sql(schema))
            if (
                schema.name not in existing_tables
                or schema_definition in self._ensured_schema_definitions
            ):
                continue  # table doesn't exist yet, or was already checked
            self._validate_table(self.connection, schema)
            ensured_definitions.append(schema_definition)
        self._ensured_schema_definitions.update(ensured_definitions)
        return existing_tables

    def _validate_table(self, connection, schema):
        """Require an existing table to match its current ixdat description."""
        table_info = {
            row[1]: {"type": row[2].upper(), "not_null": row[3], "pk": row[5]}
            for row in connection.execute(
                f"PRAGMA table_info({_quote_identifier(schema.name)})"
            ).fetchall()
        }
        if isinstance(schema, relational.LinkerTableSchema):
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
                self._incompatible_schema(
                    schema.name,
                    f"expected column {column.name!r} is missing",
                )
            actual = table_info[column.name]
            expected_type = SQLITE_TYPES.get(column.dtype, "")
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
        if column.dtype is np.ndarray:
            array = np.asarray(value)
            # Pandas represents text columns as object arrays. NumPy can save the
            # same strings without pickle once they use its native Unicode type.
            if array.dtype == object and all(
                isinstance(item, str) for item in array.flat
            ):
                array = array.astype(str)
            buffer = BytesIO()
            np.save(buffer, array, allow_pickle=False)
            return sqlite3.Binary(buffer.getvalue())
        if column.dtype in (dict, list, tuple):
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
                f"arrays, give it the type dict/list or numpy.ndarray in "
                "the `column_types` of the class which defines the column."
            )
        return value

    def _decode(self, column, value):
        """Return the python representation of a value, per its column's dtype"""
        if value is None:
            return None
        if column.dtype is np.ndarray:
            return np.load(BytesIO(value), allow_pickle=False)
        if column.dtype in (dict, list, tuple):
            decoded = json.loads(value)
            return tuple(decoded) if column.dtype is tuple else decoded
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
        definition += " " + SQLITE_TYPES[column.dtype]
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
