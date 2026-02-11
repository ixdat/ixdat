# Demo for reading and plotting a dataset from Asimov.

import argparse

import matplotlib.pyplot as plt
import numpy as np

from ixdat import Measurement


def parse_args():
    parser = argparse.ArgumentParser(
        description="Read and plot an Asimov dataset using ixdat."
    )
    parser.add_argument(
        "--dataset-id",
        required=True,
        help="Asimov dataset UUID",
    )
    parser.add_argument(
        "--version",
        type=int,
        default=None,
        help="Explicit dataset version number to fetch (defaults to latest)",
    )
    parser.add_argument(
        "--version-id",
        default=None,
        help="Explicit dataset-version UUID to fetch",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    measurement = Measurement.read(
        args.dataset_id,
        reader="asimov",
        version=args.version,
        version_id=args.version_id,
        force_login=True,
    )

    print(f"Loaded measurement: {measurement.name}")
    print(f"Technique: {measurement.technique}")
    print(f"Series available: {sorted(measurement.series_names)}")

    asimov_meta = measurement.metadata.get("asimov", {})
    print(
        "Asimov source: "
        f"dataset_id={asimov_meta.get('dataset_id')}, "
        f"version={asimov_meta.get('dataset_version')}, "
        f"version_id={asimov_meta.get('dataset_version_id')}"
    )

    u_name = None
    for candidate in ["Ewe/V", "<Ewe>/V", "raw_potential"]:
        try:
            t_u, u = measurement.grab(candidate)
            u_name = candidate
            break
        except Exception:
            continue
    if u_name is None:
        raise KeyError(
            "Could not find potential series (tried Ewe/V, <Ewe>/V, raw_potential)."
        )

    j_name = None
    for candidate in ["<I>/mA", "I/mA", "raw_current"]:
        try:
            t_j, j = measurement.grab(candidate)
            j_name = candidate
            break
        except Exception:
            continue
    if j_name is None:
        raise KeyError(
            "Could not find current series (tried <I>/mA, I/mA, raw_current)."
        )

    if t_u.shape == t_j.shape and np.allclose(t_u, t_j):
        j_for_u = j
    else:
        j_for_u = np.interp(t_u, t_j, j)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(u, j_for_u, ".", color="C0", ms=2)
    ax.set_title(f"{measurement.name}: {j_name} vs {u_name}")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(u_name)
    ax.set_ylabel(j_name)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
