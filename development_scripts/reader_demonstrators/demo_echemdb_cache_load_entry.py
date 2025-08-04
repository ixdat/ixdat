import requests
import zipfile
import io
import json
import pandas as pd

from ixdat.tools import get_default_cache_dir

# NOTE: Requires get_default_cache_dir() from ixdat.tools (in this development branch, 'echemdb')


def load_echemdb_dataset(dataset_name, version=None):
    """
    Load all CSV + JSON entries and citation for a given dataset from echemdb.

    Args:
        dataset_name (str): e.g. "alves_2011_electrochemistry_6010"
        version (str, optional): If not provided, fetch latest from PyPI.

    Returns:
        dict: {
            'version': str,
            'entries': {
                'filename.csv': {'data': DataFrame, 'info': dict}
            },
            'bib': str or None
        }
    """
    # Get latest version db from PyPI
    if version is None:
        r = requests.get("https://pypi.org/pypi/echemdb-ecdata/json")
        r.raise_for_status()
        version = r.json()["info"]["version"]
        print(f"Using latest version: {version}")

    zip_url = f"https://github.com/echemdb/electrochemistry-data/releases/download/{version}/data-{version}.zip"
    cache_root = get_default_cache_dir("echemdb", subdir=f"data-{version}")
    prefix = f"data/generated/svgdigitizer/{dataset_name}/"

    if not (cache_root / prefix).exists():
        print("Downloading and extracting dataset...")
        response = requests.get(zip_url)
        response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            selected = [f for f in z.namelist() if f.startswith(prefix)]
            if not selected:
                raise ValueError("Dataset not found in ZIP archive.")
            for f in selected:
                z.extract(f, path=cache_root)
        print("Files extracted.")
    else:
        print("Using cached files.")

    data_dir = cache_root / prefix
    entries = {}
    for csv_file in data_dir.glob("*.csv"):
        json_file = csv_file.with_suffix(".json")
        if json_file.exists():
            with open(json_file, "r") as f:
                info = json.load(f)
        else:
            info = {}
        df = pd.read_csv(csv_file)
        entries[csv_file.name] = {"data": df, "info": info}

    bib_file = data_dir / f"{dataset_name}.bib"
    bib_text = bib_file.read_text() if bib_file.exists() else None

    return {"version": version, "entries": entries, "bib": bib_text}


if __name__ == "__main__":
    result = load_echemdb_dataset("alves_2011_electrochemistry_6010")

    print("Version:", result["version"])
    print("Entries:", list(result["entries"].keys()))
    print(
        "First entry metadata:",
        result["entries"]["alves_2011_electrochemistry_6010_f1a_solid.csv"]["info"],
    )
    print("Citation (bib):", result["bib"])
