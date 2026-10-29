"""Run the main() of every demo in ixdat.demos and report which ones fail.

Usage:
    python development_scripts/run_all_demos.py          # no plots
    python development_scripts/run_all_demos.py --show   # show each demo's plots

The demos run in a temporary folder, so the files they export are removed afterwards.
Most demos need the demo data: ``python -m ixdat.demos.download``.
"""

import importlib
import os
import pkgutil
import sys
import tempfile
import traceback

from matplotlib import pyplot as plt

import ixdat.demos


def run_all_demos(show=False):
    """Run every demo's main(show=show). Return a dict of {demo name: error}."""
    failures = {}
    n_demos = 0
    for module_info in pkgutil.walk_packages(ixdat.demos.__path__, "ixdat.demos."):
        module = importlib.import_module(module_info.name)
        if not hasattr(module, "main") or not hasattr(module, "EXAMPLES"):
            continue
        n_demos += 1
        name = module_info.name.split(".", 2)[-1]
        print(f"running {name} ...", flush=True)
        try:
            module.main(show=show)
        except Exception:
            failures[name] = traceback.format_exc(limit=1).strip().splitlines()[-1]
        plt.close("all")

    print(f"\n{n_demos - len(failures)} of {n_demos} demos ran without errors.")
    for name, error in failures.items():
        print(f"  FAILED {name}: {error}")
    return failures


if __name__ == "__main__":
    show = "--show" in sys.argv
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        failures = run_all_demos(show=show)
    sys.exit(1 if failures else 0)
