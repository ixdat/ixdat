"""Demos of ixdat, grouped like ixdat itself into readers, calculators and techniques.

Each demo module has loader functions which return the demo's ixdat objects, an
``EXAMPLES`` dictionary of these loaders, and a ``main()`` which plots the results.
Most demos read data from the folder given by :func:`get_demo_data_dir`. Download
it with ``python -m ixdat.demos.download`` (requires ``pip install pooch``).
"""

from .demo_tools import (  # noqa
    download_demo_data,
    get_demo_data_dir,
    get_test_data_dir,
    view_tables,
)
