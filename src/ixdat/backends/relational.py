"""This module derives a relational database schema from ixdat's Saveable classes

Every Saveable class already describes its own persistence declaratively:

- ``table_name`` names the main table of the class,
- ``column_attrs`` names the attributes stored as columns of the main table,
- ``extra_column_attrs`` names one-to-one *extension tables* which add columns
  for inheriting classes without changing the main table (e.g. the
  "ec_measurements" table adds "ec_technique" to rows of "measurement"), and
- ``extra_linkers`` names many-to-many *linker tables* which relate rows of the
  main table to rows of another table in an ordered way (e.g. the
  "measurement_series" table relates a measurement to its data series).

This module turns that metadata into explicit, dialect-free table descriptions
(`TableSchema` and `LinkerTableSchema`), which a relational backend (see
:class:`~ixdat.backends.sqlite_backend.SQLiteBackend`) translates into DDL and
queries. It is the successor to the table definitions proposed in
https://github.com/ixdat/ixdat/pull/75, but derives the schema from the class
attributes which ixdat already maintains instead of introducing new ones.

Class inheritance maps to the database as follows: all classes sharing a
``table_name`` (e.g. all Measurement subclasses) share the main table, and each
subclass stores its additional attributes in its extension tables, joined on the
main table's primary key. A row's concrete class is not stored explicitly; it is
recovered by the classes' own ``from_dict`` dispatch (on "technique",
"series_type" or "calculator_type"), exactly as for the directory backend.

Columns have a *logical type*, used by backends to decide how a value is encoded:

- "INTEGER", "REAL", "TEXT": scalars, stored natively,
- "JSON": dicts and lists, stored as JSON text,
- "NDARRAY": numpy arrays, stored as binary blobs and loaded lazily,
- None: unspecified; stored as whatever scalar type the value has.

The logical types of ixdat's known column attributes are collected here in
``COLUMN_TYPES``; attributes not listed are stored dynamically, so a class with
new str/int/float attributes works without any registration.
"""

COLUMN_TYPES = {
    # numpy data, stored as a blob and loaded lazily:
    "data": "NDARRAY",
    # json-serializable dicts and lists:
    "metadata": "JSON",
    "aliases": "JSON",
    "tspan_bg": "JSON",
    # known floats:
    "tstamp": "REAL",
    "bg": "REAL",
    "F": "REAL",
    "RE_vs_RHE": "REAL",
    "A_el": "REAL",
    "R_Ohm": "REAL",
    # known strings:
    "name": "TEXT",
    "technique": "TEXT",
    "ec_technique": "TEXT",
    "unit_name": "TEXT",
    "series_type": "TEXT",
    "sample_name": "TEXT",
    "calculator_type": "TEXT",
    "mol": "TEXT",
    "mass": "TEXT",
    "cal_type": "TEXT",
    # single references to rows of other tables:
    "field_id": "INTEGER",
    "spectrum_id": "INTEGER",
}

FOREIGN_KEY_COLUMNS = {
    # {column attribute: the (table, column) its value refers to}
    "field_id": ("data_series", "id"),
    "spectrum_id": ("spectrums", "id"),
}


class ColumnSchema:
    """The description of one column of a table"""

    def __init__(self, name, dtype=None, foreign_key=None):
        """Initiate a column description

        Args:
            name (str): The name of the column, which is also the name of the
                attribute of the Saveable class that it stores
            dtype (str or None): The logical type of the column (see module
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

    def __init__(self, name, columns, extends=None):
        """Initiate a table description

        Args:
            name (str): The name of the table
            columns (list of ColumnSchema): The columns, starting with "id"
            extends (str or None): The name of the main table this table
                extends, if it is an extension table
        """
        self.name = name
        self.columns = columns
        self.extends = extends

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
    """The description of a many-to-many linker table

    A linker table relates rows of an owner table to rows of a linked table,
    preserving order. For example, "measurement_series" has the columns
    (measurement_id, position, data_series_id), and each of its rows represents
    a measurement's ownership of one data series. Its rows are built from and
    loaded into the id-list attribute (e.g. "s_ids") of the owner class.
    """

    def __init__(self, name, owner_table, linked_table, id_attr):
        """Initiate a linker table description

        Args:
            name (str): The name of the linker table
            owner_table (str): The name of the owning main table
            linked_table (str): The name of the table linked to
            id_attr (str): The attribute of the owner class with the id or
                ordered id's of the linked rows. By ixdat convention, id-list
                attributes end in "_ids" (e.g. "s_ids") while single-reference
                attributes end in "_id" (e.g. "ref_id").
        """
        self.name = name
        self.owner_table = owner_table
        self.linked_table = linked_table
        self.id_attr = id_attr

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

    @property
    def is_list(self):
        """Whether the attribute holds a list of id's rather than a single id"""
        return self.id_attr.endswith("_ids")

    def __repr__(self):
        return (
            f"LinkerTableSchema('{self.name}', "
            f"{self.owner_table} -> {self.linked_table})"
        )


def _ordered(attrs):
    """Return the attribute names in deterministic order, with "name" first"""
    return sorted(attrs, key=lambda attr: (attr != "name", attr))


def _column(attr):
    """Return the ColumnSchema for a column attribute"""
    return ColumnSchema(
        attr,
        dtype=COLUMN_TYPES.get(attr),
        foreign_key=FOREIGN_KEY_COLUMNS.get(attr),
    )


def main_table_schema(cls):
    """Return the TableSchema of the main table of a Saveable class"""
    columns = [ColumnSchema("id", dtype="INTEGER")]
    columns += [_column(attr) for attr in _ordered(cls.column_attrs or [])]
    return TableSchema(cls.table_name, columns)


def extension_table_schemas(cls):
    """Return the list of TableSchema of the extension tables of a Saveable class"""
    schemas = []
    for table_name, attrs in (cls.extra_column_attrs or {}).items():
        columns = [
            ColumnSchema("id", dtype="INTEGER", foreign_key=(cls.table_name, "id"))
        ]
        columns += [_column(attr) for attr in _ordered(attrs)]
        schemas.append(TableSchema(table_name, columns, extends=cls.table_name))
    return schemas


def linker_table_schemas(cls):
    """Return the list of LinkerTableSchema of the linker tables of a Saveable class"""
    return [
        LinkerTableSchema(table_name, cls.table_name, linked_table, id_attr)
        for table_name, (linked_table, id_attr) in (cls.extra_linkers or {}).items()
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

    When loading a row, the concrete class is not known in advance (e.g.
    `Measurement.get(i)` may need to build an ECMSMeasurement), so a backend
    collects the row's attributes from every extension and linker table that any
    class sharing the main table defines. Each row only has entries in the
    tables its concrete class wrote, so no wrong attributes are picked up.

    Only classes that have been imported are found. Importing ixdat imports all
    of its own techniques, so this is only a consideration for external plugins.

    Returns:
        (list of TableSchema, list of LinkerTableSchema): The extension table
            and linker table descriptions, deduplicated by table name.
    """
    extension_schemas = {}
    linker_schemas = {}
    for family_cls in saveable_classes():
        if family_cls.table_name != cls.table_name:
            continue
        for schema in extension_table_schemas(family_cls):
            extension_schemas.setdefault(schema.name, schema)
        for schema in linker_table_schemas(family_cls):
            linker_schemas.setdefault(schema.name, schema)
    return list(extension_schemas.values()), list(linker_schemas.values())
