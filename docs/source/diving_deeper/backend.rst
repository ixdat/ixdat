.. _backend:

Database backends
=================

An ixdat object remembers the backend from which it was loaded. This makes linked
objects and numerical arrays lazy: metadata is loaded immediately, while referenced
objects and array data are read only when accessed.

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
accessed. Use ``db_path=":memory:"`` for a temporary, in-memory SQLite database.

An object can also be loaded from a backend without permanently changing the active
database::

    from ixdat.backends import SQLiteBackend

    archive = SQLiteBackend("archive.sqlite")
    try:
        old_measurement = Measurement.get(12, backend=archive)
        old_data = old_measurement.series_list[0].data
    finally:
        archive.close()

Access lazy children and arrays before closing their backend. A context manager is
convenient when the lifetime is limited::

    from ixdat.db import DB

    with SQLiteBackend("archive.sqlite") as archive:
        with DB.temporary_backend(archive):
            measurement = Measurement.load("my measurement")
            data = measurement.series_list[0].data

Persistence guarantees
----------------------

Saving a ``Saveable`` and its complete graph of child objects is one transaction. If
any row, extension table, linker, or array fails to save, none of that graph's changes
are committed and newly assigned in-memory identities are restored.

Foreign-key checks are enabled on every backend connection. Ordered relationships are
stored in linker tables, and subclass-only values are stored in one-to-one extension
tables. Frequently used name, discriminator, and reverse-link columns are indexed.

The database contains an internal ``_ixdat_metadata`` table with its storage schema
version. On first use, ixdat validates each table against the schema derived from its
``Saveable`` class. Newly introduced nullable columns are added automatically. An
incompatible column type, primary key, relationship, or newer schema version raises a
``DataBaseError`` instead of risking a misread or partial write.

Connections and threads
-----------------------

Call ``backend.close()`` when finished. Switching the global database does not close
the previous backend because objects loaded from it may still need its connection for
lazy data. An ``SQLiteBackend`` connection belongs to the thread which created it;
create a separate backend instance per thread.

Known limitation: multiply-inheriting classes and extension tables
--------------------------------------------------------------------

A ``Saveable`` class which inherits from more than one other table-defining class
(for example ``ECMSMeasurement``, which inherits from both ``ECMeasurement`` and
``MSMeasurement``) declares its own ``extra_column_attrs``, which *replaces* rather
than merges with the ``extra_column_attrs`` of each parent. In practice this means
such a class writes only to the one extension table it declares itself (e.g.
``ecms_measurements``), not to the extension tables of its other parents (e.g.
``ec_measurements``) — even though the row is, semantically, also an EC measurement.

Round-tripping is unaffected: every attribute is still saved and reloaded correctly,
because classes with this shape (``ECMSMeasurement``, ``MSSpectroMeasurement``, ...)
already duplicate the small number of overlapping columns into their own extension
table by hand. What breaks is a query that joins against only *one* parent's
extension table expecting it to be exhaustive, for example::

    -- misses every ECMSMeasurement row, though they are EC measurements too:
    SELECT m.id, m.name, e.ec_technique
    FROM measurement m
    JOIN ec_measurements e ON e.id = m.id

Filter or group by ``measurement.technique`` instead of joining a single extension
table when you want "every measurement of a given family", or ``UNION`` across the
extension tables of the classes you care about.

This is not new to the SQLite backend: it is a limitation of how ``Saveable``
resolves ``extra_column_attrs`` across multiple inheritance
(``ixdat.db.Saveable.get_main_dict``/``as_dict``), so it equally affects
``obj.as_dict()`` on the directory backend for any class which does not manually
duplicate its overlapping columns the way ``ECMSMeasurement`` does today. Properly
merging table definitions across multiple inheritance was the central open question
of `PR #75 <https://github.com/ixdat/ixdat/pull/75>`_ and remains unresolved;
fixing it belongs in ``Saveable`` itself; it is out of scope for the database
backends.

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
