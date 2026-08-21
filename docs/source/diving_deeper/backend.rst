.. _backend:

Database backends
=================

An ixdat object remembers the backend from which it was loaded. This keeps linked
objects and ``DataSeries.data`` out of the first load. ixdat reads them from the
remembered backend when they are first accessed.

SQLite
------

The SQLite backend stores a complete ixdat project in one portable ``.sqlite`` file.
It uses Python's built-in :mod:`sqlite3` module and requires no optional dependency::

    from ixdat.db import change_database

    backend = change_database("sqlite", db_path="my_project.sqlite")
    measurement_id = measurement.save()

    by_id = Measurement.get(measurement_id)
    by_name = Measurement.load(measurement.name)

``Saveable.load(name)`` returns the most recently saved row with that exact name.
``Saveable.get(id)`` selects a specific row. Arrays stored in data series are NumPy
``.npy`` payloads inside SQLite blob columns and remain unloaded until ``.data`` is
accessed. This delayed loading belongs specifically to the column named ``data``.
Another NumPy-array column, such as one added by a plugin, loads with the object's
other stored values.

The string ``":memory:"`` is a special SQLite address. Calling
``sqlite3.connect(":memory:")`` creates a temporary database which belongs to that
one connection and has no file on disk. Calling it again creates another empty
database, even though both calls use the same address string. ixdat therefore treats
each in-memory backend as a separate database and gives it a unique internal address.
This keeps objects from two temporary databases from receiving the same identity when
both happen to contain, for example, a ``data_series`` row with ``id=1``. Closing the
connection removes its temporary database. Use ``db_path=":memory:"`` when this
short-lived storage is useful.

Two ``SQLiteBackend`` objects can open separate connections to the same file. They
remain separate Python objects, and ``shares_storage_with()`` reports that both
connections reach the same saved rows. ixdat uses this explicit check when deciding
whether an object already belongs to the target database file. SQLite coordinates
writes made through the separate connections.

An object's ``id`` is its integer row number inside one backend. Its
``short_identity`` is always ``(backend, id)``, so a reference still identifies the
right object when one measurement uses data from several backends. A storage backend
turns these references into local integer ids when it saves the connected objects.
The backend in this tuple is the live Python backend object. Two backend objects can
open separate connections to the same file, so compare such references with
``same_short_identity()`` when their table is already known::

    from ixdat.db import same_short_identity

    same_short_identity(first.short_identity, second.short_identity)

An object can also be loaded from a backend without permanently changing the active
database::

    from ixdat.backends import SQLiteBackend

    archive = SQLiteBackend("archive.sqlite")
    try:
        old_measurement = Measurement.get(12, backend=archive)
        old_data = old_measurement.series_list[0].data
    finally:
        archive.close()

Every backend can be used as a context manager. Leaving the ``with`` block calls
``backend.close()``. For SQLite, this closes the open connection. Access lazily
loaded related objects and ``.data`` arrays before leaving the block::

    from ixdat.db import DB

    with SQLiteBackend("archive.sqlite") as archive:
        with DB.temporary_backend(archive):
            measurement = Measurement.load("my measurement")
            data = measurement.series_list[0].data

Persistence guarantees
----------------------

Saving a ``Saveable`` and all the related objects it owns is one transaction. If any
row, extension table, table of connections, or array fails to save, none of those
changes are committed and newly assigned in-memory identities are restored.

Foreign-key checks are enabled on every backend connection. Ordered relationships are
stored in tables of connections, also called linker tables. Values added by a more
specialized class are stored in one-to-one extension tables. Frequently used name,
discriminator, and reverse-link columns are indexed.

Every new database contains an internal ``_ixdat_metadata`` table recording storage
schema version 1. ixdat opens an existing database only when it records the same
version. An unversioned database or a different version raises ``DataBaseError`` with
instructions to migrate its objects into a new database.

Saving an object creates any tables needed for object types that have not previously
been saved. Reading validates existing tables and leaves the database unchanged. A
missing column, incompatible column type, primary key, or relationship raises
``DataBaseError`` and requires an explicit migration. A future schema upgrade should
read objects from a supported old database and save them into a new database, leaving
the source file unchanged.

Defining the tables of a ``Saveable`` class
-------------------------------------------

The relational schema is not written by hand. It is derived, by
``ixdat.backends.relational``, from the table metadata every ``Saveable`` class
already carries. ``table_name`` and ``column_attrs`` define ordinary values in the
main table, and ``extra_column_attrs`` defines ordinary values in extension tables.
Relationships to other saved objects are declared with ``Relationship``. These
definitions are merged over class ancestry, so each class only declares what it adds.

For example, this class stores two ordinary values and one reference to a saved
calculator::

    from ixdat.db import Relationship

    class XRFSpectrum(Spectrum):
        extra_column_attrs = {
            "xrf_spectrums": {"excitation_energy", "detector_gains"}
        }
        column_types = {"excitation_energy": float, "detector_gains": list}
        relationships = {
            "calibration": Relationship(
                "calculator",
                "calibration_id",
                storage_table="xrf_spectrums",
            )
        }

``column_types`` gives the Python type of a column's value. Supported declarations
are ``int``, ``float``, ``str``, ``dict``, ``list``, ``tuple``, and
``numpy.ndarray``. Each relational backend maps these to its own storage types.
SQLite stores containers as JSON and arrays as blobs. As described above, the column
named ``data`` stays in the file until its property is accessed; other array columns
load with the object's ordinary values. A column which is not listed has no fixed
type and stores its scalar value directly. Saving a container or array to an untyped
column raises a ``DataBaseError`` naming the column.

The key in ``relationships`` is the Python attribute holding the linked object.
``Relationship`` names the linked table and the attribute holding its saved id.
``storage_table`` says where ixdat stores that id. Leaving it out for a single object
places the id in the main table. A list uses ``many=True`` and requires a
``storage_table`` for its rows of connections. Each row receives a ``position`` so
ixdat can rebuild the list in the same order, including repeated objects.

By default, saving a new owner saves new related objects first. This gives each related
object the id which the owner needs to store. A normal repeated ``save()`` keeps
already-saved rows unchanged. After changing a saved object or replacing one of its
relationships, ``DB.backend.save(owner, force=True)`` updates the owner and its related
objects together. Set ``save_related=False`` when the related object should stay
outside the owner's save process.

Older plugin classes can continue to use ``extra_linkers``, ``column_references``, and
``child_attrs``. ixdat converts those declarations into the same table descriptions.
For an older ``extra_linkers`` declaration, an id attribute ending in ``_ids`` means
that it contains an ordered list.

Keeping this information on the classes lets a ``Saveable`` class in a plugin or
user script define its own relational tables directly.

Extension tables represent inheritance without repeating parent columns. For example,
``ECMSMeasurement`` contributes ``tspan_bg`` to ``ecms_measurements`` and inherits
``ec_technique`` from ``ECMeasurement`` in ``ec_measurements``. Both tables use the
measurement's main-table id as their primary and foreign key.

Connections and threads
-----------------------

All backends provide ``close()`` and support the context-manager form shown above. A
backend overrides ``close()`` when it owns a long-lived resource. SQLite uses it to
close the open connection. Switching the global database does not close the previous
backend because objects loaded from it may still need its connection for lazy data.
An ``SQLiteBackend`` connection belongs to the thread which created it; create a
separate backend instance per thread.

The directory backend
---------------------

The directory backend stores JSON metadata and NumPy arrays in a folder hierarchy::

    change_database(
        "directory",
        directory=Path("data"),
        project_name="my_project",
    )

It implements the same ``save()``, ``get(id)``, and exact ``load(name)`` interface,
but SQLite is preferable when transactional saves, relational queries, compact
distribution, and schema validation are important.
