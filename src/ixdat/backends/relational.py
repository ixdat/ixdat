"""Derive backend-neutral table descriptions from Saveable metadata.

This module converts the persistence information on Saveable classes into
ColumnSchema, TableSchema, and LinkerTableSchema objects for relational
backends.

See :ref:`backend` for the persistence model, supported types, inheritance,
and linked objects.
"""

import numpy as np

from ..exceptions import DataBaseError


KNOWN_COLUMN_TYPES = (int, float, str, dict, list, tuple, np.ndarray)


def _type_name(dtype):
    """Return a readable name for a declared Python column type."""
    if isinstance(dtype, str):
        return dtype
    return getattr(dtype, "__name__", repr(dtype))


class ColumnSchema:
    """The description of one column of a table"""

    def __init__(self, name, dtype=None, foreign_key=None):
        """Initiate a column description

        Args:
            name (str): The name of the column, which is also the name of the
                attribute of the Saveable class that it stores
            dtype (type or None): The Python type of the column (see module
                docstring). None means unspecified/dynamic.
            foreign_key (tuple or None): The (table_name, column_name) that this
                column's values refer to, if it is a foreign key
        """
        self.name = name
        self.dtype = dtype
        self.foreign_key = foreign_key

    def __repr__(self):
        return f"ColumnSchema('{self.name}', dtype={self.dtype})"

    def __eq__(self, other):
        return (
            self.__class__ is other.__class__
            and self.name == other.name
            and self.dtype == other.dtype
            and self.foreign_key == other.foreign_key
        )


class TableSchema:
    """The description of a main table or a one-to-one extension table

    A main table has one row per object and an "id" primary key column. An
    extension table also has one row per object, but its "id" column is both
    primary key and foreign key to the main table it extends - this is how
    columns of inheriting classes are represented without changing the main
    table (see module docstring).
    """

    def __init__(self, name, columns, base_table=None):
        """Initiate a table description

        Args:
            name (str): The name of the table
            columns (list of ColumnSchema): The columns, starting with "id"
            base_table (str or None): The name of the main table this extension
                table adds values to
        """
        self.name = name
        self.columns = columns
        self.base_table = base_table

    @property
    def column_names(self):
        return [column.name for column in self.columns]

    @property
    def data_columns(self):
        """The columns which hold attribute values, i.e. all but "id" """
        return [column for column in self.columns if column.name != "id"]

    def __repr__(self):
        return f"TableSchema('{self.name}', columns={self.column_names})"


class LinkerTableSchema:
    """The description of a table which stores connections between saved objects

    A linker table relates rows of an owner table to rows of a linked table,
    preserving order for a list. For example, "measurement_series" has the columns
    (measurement_id, position, data_series_id), and each of its rows represents
    a measurement's ownership of one data series. Its rows are built from and
    loaded into the id attribute (e.g. "s_ids") of the owner class.
    """

    def __init__(self, name, owner_table, linked_table, id_attr, many):
        """Initiate a linker table description

        Args:
            name (str): The name of the linker table
            owner_table (str): The name of the owning main table
            linked_table (str): The name of the table linked to
            id_attr (str): The attribute of the owner class with the id or
                ordered ids of the linked rows
            many (bool): Whether the id attribute holds an ordered list. New
                ``Relationship`` definitions provide this directly. Older
                ``extra_linkers`` definitions use their established ``_ids`` naming
                rule.
        """
        self.name = name
        self.owner_table = owner_table
        self.linked_table = linked_table
        self.id_attr = id_attr
        self.many = many

    @property
    def owner_column(self):
        """The name of the column with the id of the owning row"""
        return self.owner_table + "_id"

    @property
    def linked_column(self):
        """The name of the column with the id of the linked row"""
        if self.linked_table == self.owner_table:
            # e.g. "component_measurements", which links measurements to
            # measurements, gets the columns (measurement_id, linked_measurement_id)
            return "linked_" + self.linked_table + "_id"
        return self.linked_table + "_id"

    def __repr__(self):
        return (
            f"LinkerTableSchema('{self.name}', "
            f"{self.owner_table} -> {self.linked_table})"
        )


def _ordered(attrs):
    """Return the attribute names in deterministic order, with "name" first"""
    return sorted(attrs, key=lambda attr: (attr != "name", attr))


def _validate_column_metadata(cls):
    """Reject type or reference declarations for columns the class does not store."""
    stored_attrs = set(cls.get_main_column_attrs())
    for attrs in cls.get_extra_column_attrs().values():
        stored_attrs.update(attrs)

    declarations = {
        "column_types": cls.get_column_types(),
        "column_references": cls.get_column_references(),
    }
    for metadata_name, metadata in declarations.items():
        unknown_attrs = set(metadata) - stored_attrs
        if unknown_attrs:
            raise DataBaseError(
                f"{cls.__name__} declares {metadata_name} for unknown column(s): "
                f"{', '.join(sorted(unknown_attrs))}."
            )
    for attr, dtype in declarations["column_types"].items():
        if dtype is not None and dtype not in KNOWN_COLUMN_TYPES:
            raise DataBaseError(
                f"{cls.__name__} gives its column '{attr}' the unknown type "
                f"{_type_name(dtype)!r}. The types a column can have are "
                f"{', '.join(_type_name(known) for known in KNOWN_COLUMN_TYPES)}."
            )


def _columns(cls, attrs):
    """Return the ColumnSchemas for the named column attributes of a Saveable class"""
    types = cls.get_column_types()
    references = cls.get_column_references()
    columns = []
    for attr in _ordered(attrs):
        dtype = types.get(attr)
        referenced_table = references.get(attr)
        columns.append(
            ColumnSchema(
                attr,
                # a column holding another row's id is that id's type:
                dtype=int if referenced_table else dtype,
                foreign_key=(referenced_table, "id") if referenced_table else None,
            )
        )
    return columns


def main_table_schema(cls):
    """Return the TableSchema of the main table of a Saveable class"""
    _validate_column_metadata(cls)
    columns = [ColumnSchema("id", dtype=int)]
    columns += _columns(cls, cls.get_main_column_attrs())
    return TableSchema(cls.table_name, columns)


def extension_table_schemas(cls):
    """Return the extension-table schemas of a Saveable class and its ancestors."""
    _validate_column_metadata(cls)
    schemas = []
    for table_name, attrs in cls.get_extra_column_attrs().items():
        columns = [ColumnSchema("id", dtype=int, foreign_key=(cls.table_name, "id"))]
        columns += _columns(cls, attrs)
        schemas.append(TableSchema(table_name, columns, base_table=cls.table_name))
    return schemas


def linker_table_schemas(cls):
    """Return the connection-table schemas of a Saveable class and its ancestors.

    New relationships state ``many`` directly. Older ``extra_linkers`` declarations
    keep their established rule: an id attribute ending in ``_ids`` contains a list.
    """
    many_by_table = {
        relationship.storage_table: relationship.many
        for relationship in cls.get_relationships().values()
        if relationship.many
    }
    return [
        LinkerTableSchema(
            table_name,
            cls.table_name,
            linked_table,
            id_attr,
            many=many_by_table.get(table_name, id_attr.endswith("_ids")),
        )
        for table_name, (linked_table, id_attr) in cls.get_extra_linkers().items()
    ]


def table_schemas_of(cls):
    """Return all table descriptions needed to save objects of a Saveable class"""
    return (
        [main_table_schema(cls)]
        + extension_table_schemas(cls)
        + linker_table_schemas(cls)
    )


def saveable_classes():
    """Return the list of all imported Saveable classes"""
    from ..db import Saveable  # here to avoid circular import

    classes = []

    def collect(cls):
        for subclass in cls.__subclasses__():
            if subclass not in classes:
                classes.append(subclass)
            collect(subclass)

    collect(Saveable)
    return classes


def family_table_schemas(cls):
    """Return the extension and linker tables of all classes sharing cls's main table

    Before reading a row, its exact subclass is unknown (a row in "measurement"
    could be an ECMSMeasurement, for example). A backend therefore checks the
    extension/linker tables of every class that shares this main table. A row has
    entries only in the tables its concrete class wrote to.

    Only classes that have already been imported can be found this way. Importing
    ixdat imports all of its own techniques, so this only matters for external
    plugin classes that haven't been imported yet.

    A class whose own table metadata is broken is skipped here, so that one
    misdeclared class doesn't stop every other row of the shared main table from
    being read. Its error is raised when that class is itself saved or loaded.

    Returns:
        (list of TableSchema, list of LinkerTableSchema): The extension table
            and linker table descriptions, deduplicated by table name.
    """
    extension_schemas = {}
    linker_schemas = {}
    for family_cls in saveable_classes():
        if family_cls.table_name != cls.table_name:
            continue
        try:
            found_extensions = extension_table_schemas(family_cls)
            found_linkers = linker_table_schemas(family_cls)
        except DataBaseError:
            if family_cls is cls:
                raise  # cls is the class we were asked about, so it must be right
            continue  # a sibling of cls; the row being read doesn't need it
        for schema in found_extensions:
            extension_schemas.setdefault(schema.name, schema)
        for schema in found_linkers:
            linker_schemas.setdefault(schema.name, schema)
    return list(extension_schemas.values()), list(linker_schemas.values())
