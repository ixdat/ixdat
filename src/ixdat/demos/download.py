"""Download the demo data: ``python -m ixdat.demos.download [cache_dir]``

Requires the `pooch` package: ``pip install pooch``.
"""

import sys

from ixdat.demos.demo_tools import download_demo_data

if __name__ == "__main__":
    data_dir = download_demo_data(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"The demo data is in {data_dir}")
