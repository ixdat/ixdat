import json, yaml
from pathlib import Path

data_dir = Path(__file__).parent.parent / "src/ixdat/plugin_data/ms_quant"

for folder_name in ["molecules", "chips"]:
    for file in (data_dir / folder_name).iterdir():
        with open(file) as f:
            data = json.load(f)
        with open(file.with_suffix(".yml"), "w") as g:
            yaml.dump(data, g)
