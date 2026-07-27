"""This module derives a relational database schema from ixdat's Saveable classes

Every Saveable class already says how it should be saved:

- ``table_name`` is the name of its main table,
- ``column_attrs`` are the attributes stored as columns of that table,
- ``extra_column_attrs`` names *extension tables*: one extra table per subclass
  that adds its own columns without changing the shared main table (e.g. the
  "ec_measurements" table adds "ec_technique" to rows of "measurement"), and
- ``extra_linkers`` names *linker tables*, which record ordered references from
  a row to several rows of another table (e.g. "measurement_series" records
  which data series belong to which measurement, and in what order),
- ``column_types`` gives the *logical type* of the columns where it matters, and
- ``column_references`` names the columns which hold the id of a row of another
  table (e.g. a spectrum's "field_id" refers to the "data_series" table).

This module turns that information into plain table descriptions (`TableSchema`
and `LinkerTableSchema`) that don't depend on any particular database - a
backend (see :class:`~ixdat.backends.sqlite_backend.SQLiteBackend`) then turns
those into real SQL. This does the same job PR #75
(https://github.com/ixdat/ixdat/pull/75) set out to do, but by reading the
class attributes ixdat already has instead of asking every class to redeclare
its columns in a new format.

Classes that share a ``table_name`` (e.g. every Measurement subclass) share one
main table; each subclass's own extra attributes go in its own extension table,
linked back to the main table's id. Which subclass a row belongs to isn't
stored as a separate column - it's worked out the same way the directory
backend already does it, from a column like "technique" or "calculator_type".

A column's *logical type* tells a backend how to store its value:

- "INTEGER", "REAL", "TEXT": plain numbers/text, stored as-is,
- "JSON": dicts and lists, stored as JSON text,
- "NDARRAY": numpy arrays, stored as a binary blob and only loaded when needed,
- None: no fixed type; stored as whatever type the value already is.

A class states these in ``column_types``, alongside the ``column_attrs`` that
introduce the columns, and they are merged over the class's ancestry by
``Saveable.get_column_types()``. A column left out still works - it's just
stored using whatever scalar type Python gives it - so a new class with a plain
str/int/float attribute needs to declare nothing. ``column_references``, merged
the same way by ``Saveable.get_column_references()``, says which columns hold
the id of a row of another table, and thus become foreign keys.

Keeping this on the classes, rather than in a table here, is what lets a class
outside of ixdat define columns of any type without ixdat knowing about it.
"""

from ..exceptions import DataBaseError


KNOWN_COLUMN_TYPES = ("INTEGER", "REAL", "TEXT", "JSON", "NDARRAY")


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


def _columns(cls, attrs):
    """Return the ColumnSchemas for the named column attributes of a Saveable class"""
    types = cls.get_column_types()
    references = cls.get_column_references()
    columns = []
    for attr in _ordered(attrs):
        dtype = types.get(attr)
        if dtype and dtype not in KNOWN_COLUMN_TYPES:
            raise DataBaseError(
                f"{cls.__name__} gives its column '{attr}' the unknown type "
                f"'{dtype}'. The types a column can have are "
                f"{', '.join(KNOWN_COLUMN_TYPES)}."
            )
        referenced_table = references.get(attr)
        columns.append(
            ColumnSchema(
                attr,
                # a column holding another row's id is that id's type:
                dtype="INTEGER" if referenced_table else dtype,
                foreign_key=(referenced_table, "id") if referenced_table else None,
            )
        )
    return columns


def main_table_schema(cls):
    """Return the TableSchema of the main table of a Saveable class"""
    columns = [ColumnSchema("id", dtype="INTEGER")]
    columns += _columns(cls, cls.column_attrs or [])
    return TableSchema(cls.table_name, columns)


def extension_table_schemas(cls):
    """Return the list of TableSchema of the extension tables of a Saveable class

    Only looks at what `cls` itself declares in `extra_column_attrs`, not what its
    parent classes declare. Usually fine, since a class only overrides this if it
    needs its own extra columns - but a class inheriting from *two* table-defining
    classes (like `ECMSMeasurement`) ends up only writing to its own extension
    table, not its other parent's. See the "Known limitation" section of
    docs/source/diving_deeper/backend.rst.
    """
    schemas = []
    for table_name, attrs in (cls.extra_column_attrs or {}).items():
        columns = [
            ColumnSchema("id", dtype="INTEGER", foreign_key=(cls.table_name, "id"))
        ]
        columns += _columns(cls, attrs)
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

    Before reading a row, we don't yet know which subclass it actually is (a row
    in "measurement" could turn out to be an ECMSMeasurement, say), so a backend
    needs to check the extension/linker tables of every class that could share
    this main table, not just `cls`. A row will only actually have entries in the
    tables its real class wrote to, so checking extra tables is harmless.

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
