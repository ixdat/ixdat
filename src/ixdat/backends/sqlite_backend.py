"""This module implements an SQLite database backend for ixdat

The SQLite backend saves ixdat objects to, and loads them from, a proper
relational database in a single local file. The schema is not written by hand:
it is derived from the table metadata that every Saveable class already defines
(see :module:`~ixdat.backends.relational`). Tables are created on demand with
``CREATE TABLE IF NOT EXISTS``, so any Saveable class - including ones defined
in external plugins - can be saved without registration here.

Use it like any other ixdat backend::

    from ixdat.db import change_database

    change_database("sqlite", db_path="my_project.sqlite")
    measurement.save()
    ...
    loaded = Measurement.get(1)
    # or, by name:
    loaded = Measurement.load("my measurement")

Numerical data (numpy arrays) is stored in NDARRAY (numpy-format blob) columns,
which are excluded from queries until the data is actually needed - keeping
ixdat's laziness intact: `Measurement.get()` touches only metadata, and each
DataSeries fetches its array from the database the first time `.data` is
accessed.

Note that an SQLiteBackend's connection belongs to the thread that created the
backend, as is the default for python's sqlite3 module.
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
        if db_path:
            self.db_path = Path(db_path)
        else:
            directory = Path(directory or config.standard_data_directory)
            project_name = project_name or config.default_project_name
            self.db_path = directory / (project_name + DATABASE_FILE_SUFFIX)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.db_path)
        self.connection.execute("PRAGMA foreign_keys = ON")
        self._ensured_classes = set()
        super().__init__()

    @property
    def address(self):
        """The path to the SQLite database file"""
        return str(self.db_path)

    def close(self):
        """Close the connection to the database file"""
        self.connection.close()

    def __eq__(self, other):
        """Two SQLiteBackends are equivalent if they use the same database file"""
        if other is self:
            return True
        return (
            other.__class__ is self.__class__
            and other.db_path.resolve() == self.db_path.resolve()
        )

    # ------- saving  ------- #

    def save(self, obj, force=False, no_updates=True):
        """Save a Saveable object as rows of the tables defined by its class

        Args:
            obj (Saveable): the object to save
            force (bool): Whether to force updates if the object is already saved
            no_updates (bool): Whether to allow updates if the object is already
                saved. If both force and no_updates are False, the user will be
                prompted on whether to save.
        """
        # First, save the objects that this object references, so that they get
        # their id's in this backend for the object's rows to correctly refer
        # to. This is done recursively, and mirrors the directory backend.
        if obj.child_attrs:
            for child_list_name in obj.child_attrs:
                child_list = getattr(obj, child_list_name) or []
                for child_obj in child_list:
                    self.save(child_obj, force=force, no_updates=True)
        self._ensure_tables(type(obj))
        # If it's already saved here, decide whether to update it:
        if obj.backend is self and self.contains(obj.table_name, obj.id):
            okay_to_update = not no_updates
            update_the_row = force or (
                okay_to_update
                and prompt_for_permission(
                    f"Are you sure you would like to overwrite "
                    f"{self} table={obj.table_name} id={obj.id} with {obj}? "
                    f"(You can use save() with force=True to suppress this.)"
                )
            )
            if update_the_row:
                self.update_row(obj)
                return obj.id  # return the id of the updated row
            else:
                return  # return nothing since nothing was done
        i = self.add_row(obj)
        obj.set_id(i)
        obj.set_backend(self)
        return i

    def add_row(self, obj):
        """Insert the object's rows in its main, extension, and linker tables"""
        main_schema = relational.main_table_schema(type(obj))
        main_dict = obj.get_main_dict()
        columns = main_schema.data_columns
        values = [self._encode(column, main_dict[column.name]) for column in columns]
        with self.connection as connection:
            cursor = connection.execute(
                _insert_sql(main_schema.name, [column.name for column in columns]),
                values,
            )
            i = cursor.lastrowid
            self._insert_extras(connection, obj, i)
        return i

    def update_row(self, obj):
        """Update the object's rows in its main, extension, and linker tables"""
        main_schema = relational.main_table_schema(type(obj))
        main_dict = obj.get_main_dict()
        columns = main_schema.data_columns
        values = [self._encode(column, main_dict[column.name]) for column in columns]
        assignments = ", ".join(f'"{column.name}" = ?' for column in columns)
        with self.connection as connection:
            connection.execute(
                f'UPDATE "{main_schema.name}" SET {assignments} WHERE "id" = ?',
                values + [obj.id],
            )
            for schema in relational.extension_table_schemas(type(obj)):
                connection.execute(
                    f'DELETE FROM "{schema.name}" WHERE "id" = ?', (obj.id,)
                )
            for linker in relational.linker_table_schemas(type(obj)):
                connection.execute(
                    f'DELETE FROM "{linker.name}" WHERE "{linker.owner_column}" = ?',
                    (obj.id,),
                )
            self._insert_extras(connection, obj, obj.id)

    def _insert_extras(self, connection, obj, i):
        """Insert the object's rows in its extension and linker tables"""
        for schema in relational.extension_table_schemas(type(obj)):
            columns = schema.data_columns
            values = [
                self._encode(column, getattr(obj, column.name)) for column in columns
            ]
            connection.execute(
                _insert_sql(
                    schema.name, ["id"] + [column.name for column in columns]
                ),
                [i] + values,
            )
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

        The object's attributes are collected from cls's main table and from the
        extension and linker tables of every class sharing that main table (only
        the tables its concrete class wrote have rows with its id). Numerical
        data is not loaded here - it is loaded lazily via `load_obj_data`.
        """
        self._ensure_tables(cls)
        main_schema = relational.main_table_schema(cls)
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
        existing_tables = self._existing_tables()
        extension_schemas, linker_schemas = relational.family_table_schemas(cls)
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
                + f' WHERE "{linker.owner_column}" = ? ORDER BY "position"',
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
            f'SELECT "id" FROM "{cls.table_name}" WHERE "name" = ? '
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
            f'SELECT 1 FROM "{table_name}" WHERE "id" = ?', (i,)
        ).fetchone()
        return row is not None

    def get_next_available_id(self, table_name, obj=None):
        """Return the next available id for a given table"""
        if table_name not in self._existing_tables():
            return 1
        row = self.connection.execute(
            f'SELECT COALESCE(MAX("id"), 0) + 1 FROM "{table_name}"'
        ).fetchone()
        return row[0]

    # ------- schema  ------- #

    def _ensure_tables(self, cls):
        """Create the tables needed by a Saveable class, if they don't yet exist"""
        if cls in self._ensured_classes:
            return
        with self.connection as connection:
            for schema in relational.table_schemas_of(cls):
                connection.execute(_create_table_sql(schema))
        self._ensured_classes.add(cls)

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
        """
        statements = {}
        for cls in relational.saveable_classes():
            if not cls.table_name:
                continue
            for schema in relational.table_schemas_of(cls):
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
            return json.dumps(value)
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (dict, list, tuple, np.ndarray)):
            raise DataBaseError(
                f"Can't save value {value!r} in the dynamically typed column "
                f"'{column.name}'. If this column should hold dicts/lists or "
                "numpy arrays, add it to COLUMN_TYPES in ixdat.backends.relational "
                "as 'JSON' or 'NDARRAY', respectively."
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
        """Turn a `short_identity` into an id, requiring that it refers to self"""
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
            f'CREATE TABLE IF NOT EXISTS "{schema.name}" (\n'
            f'    "{schema.owner_column}" INTEGER NOT NULL '
            f'REFERENCES "{schema.owner_table}"("id"),\n'
            f'    "position" INTEGER NOT NULL,\n'
            f'    "{schema.linked_column}" INTEGER NOT NULL '
            f'REFERENCES "{schema.linked_table}"("id"),\n'
            f'    PRIMARY KEY ("{schema.owner_column}", "position")\n'
            ")"
        )
    column_definitions = []
    for column in schema.columns:
        definition = f'"{column.name}"'
        if column.dtype:
            definition += " " + column.dtype
        if column.name == "id":
            definition += " PRIMARY KEY"
        if column.foreign_key:
            foreign_table, foreign_column = column.foreign_key
            definition += f' REFERENCES "{foreign_table}"("{foreign_column}")'
        column_definitions.append("    " + definition)
    return (
        f'CREATE TABLE IF NOT EXISTS "{schema.name}" (\n'
        + ",\n".join(column_definitions)
        + "\n)"
    )


def _insert_sql(table_name, column_names):
    """Return a parametrized INSERT statement for the named table and columns"""
    columns = ", ".join(f'"{name}"' for name in column_names)
    placeholders = ", ".join("?" for _ in column_names)
    return f'INSERT INTO "{table_name}" ({columns}) VALUES ({placeholders})'


def _select_sql(table_name, column_names):
    """Return a SELECT statement for the named table and columns"""
    columns = ", ".join(f'"{name}"' for name in column_names)
    return f'SELECT {columns} FROM "{table_name}"'
