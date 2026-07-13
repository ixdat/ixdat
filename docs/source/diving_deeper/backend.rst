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
