"""Save the examples of all the reader demos in an SQLite file and view its tables.

Runs every loader in the ``EXAMPLES`` of every demo module in ``ixdat.demos``, saves the
returned objects, and prints the examples which could not be loaded or saved.
"""

import importlib
import pkgutil
import traceback
from pathlib import Path

import ixdat.demos
from ixdat.db import change_database
from ixdat.demos import view_tables

sqlite_file = Path(".") / "all_demo_data.sqlite"
if sqlite_file.exists():
    sqlite_file.unlink()

change_database("sqlite", db_path=sqlite_file)

failures = {}
for module_info in pkgutil.walk_packages(ixdat.demos.__path__, "ixdat.demos."):
    module = importlib.import_module(module_info.name)
    for name, load in getattr(module, "EXAMPLES", {}).items():
        key = f"{module_info.name.split('.', 2)[-1]}.{name}"
        try:
            loaded = load()
            for obj in loaded if isinstance(loaded, list) else [loaded]:
                obj.save()
        except Exception:
            failures[key] = traceback.format_exc(limit=1).strip().splitlines()[-1]
            continue
        print(f"saved {key}")

print(f"\n{len(failures)} examples failed:")
for key, error in failures.items():
    print(f"  {key}: {error}")

view_tables(sqlite_file)
