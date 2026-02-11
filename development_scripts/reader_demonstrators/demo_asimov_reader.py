"""Demo for reading and plotting a dataset from Asimov with Keycloak auth.

Usage:
    python development_scripts/reader_demonstrators/demo_asimov_reader.py \
        --dataset-id <DATASET_UUID>

Environment variables (recommended):
    ASIMOV_BASE_URL=https://asimov.enci.dk/api
    KEYCLOAK_CLIENT_SECRET=<optional-for-confidential-client>
    ASIMOV_ACCESS_TOKEN=<optional-bearer-token>

First run opens a browser window for device login. Tokens are cached locally
for later runs (unless --force-login is passed).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ixdat import Measurement
from ixdat.readers.asimov import AsimovReader


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
    parser.add_argument(
        "--token",
        default=None,
        help="Bearer token (optional; otherwise Keycloak flow/env is used)",
    )
    parser.add_argument(
        "--base-url",
        default=None,
        help="Asimov API base URL (default: env ASIMOV_BASE_URL or https://asimov.enci.dk/api)",
    )
    parser.add_argument(
        "--keycloak-client-secret",
        default=None,
        help="Keycloak client secret (default: env KEYCLOAK_CLIENT_SECRET)",
    )
    parser.add_argument(
        "--force-login",
        action="store_true",
        help="Force a fresh Keycloak device login",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not auto-open browser in Keycloak device flow",
    )
    parser.add_argument(
        "--savefig",
        default=None,
        help="Optional path to save figure image",
    )
    return parser.parse_args()


def plot_measurement(measurement):
    try:
        axes = measurement.plot()
        if axes:
            return axes[0].get_figure()
    except Exception as exc:
        print(f"Standard ixdat plot failed ({exc}); using manual fallback plot.")

    value_names = sorted(list(measurement.value_names))
    if not value_names:
        raise RuntimeError("Measurement has no value series for manual fallback plot.")
    series_name = value_names[0]

    t, v = measurement.grab(series_name)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, v, color="C0", lw=1.5)
    ax.set_title(measurement.name)
    ax.set_xlabel(measurement.t_name)
    ax.set_ylabel(series_name)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def plot_ewe_vs_current(measurement):
    series_names = measurement.series_names
    x_name = "Ewe/V" if "Ewe/V" in series_names else "raw_potential"
    y_name = "<I>/mA" if "<I>/mA" in series_names else "raw_current"

    try:
        t_x, x = measurement.grab(x_name)
        t_y, y_raw = measurement.grab(y_name)
    except Exception as exc:
        print(
            f"Skipping Ewe-vs-I plot: could not resolve '{x_name}'/'{y_name}' ({exc})."
        )
        return None

    if t_x.shape == t_y.shape and np.allclose(t_x, t_y):
        y = y_raw
        pairing = "direct (shared time axis)"
    else:
        # If axes differ, align current to potential time (same strategy as ECPlotter).
        y = np.interp(t_x, t_y, y_raw)
        pairing = "interpolated onto potential time axis"

    print(f"Ewe-vs-I plot uses x='{x_name}', y='{y_name}', pairing={pairing}.")

    fig, ax = plt.subplots(figsize=(6, 5))
    # Scatter emphasizes measured points and avoids misleading visual continuity.
    ax.plot(x, y, ".", color="C1", ms=2, alpha=0.8)
    ax.set_title(f"{measurement.name}: {y_name} vs {x_name}")
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def main():
    args = parse_args()

    reader = AsimovReader(
        base_url=args.base_url,
        token=args.token,
        keycloak_client_secret=args.keycloak_client_secret,
        keycloak_open_browser=(not args.no_browser),
    )

    measurement = Measurement.read(
        args.dataset_id,
        reader=reader,
        version=args.version,
        version_id=args.version_id,
        force_login=args.force_login,
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

    fig = plot_measurement(measurement)
    iv_fig = plot_ewe_vs_current(measurement)
    if args.savefig:
        save_path = Path(args.savefig)
        fig.savefig(save_path, dpi=200)
        print(f"Saved figure to: {args.savefig}")
        if iv_fig:
            iv_path = save_path.with_name(
                f"{save_path.stem}_ewe_vs_i{save_path.suffix or '.png'}"
            )
            iv_fig.savefig(iv_path, dpi=200)
            print(f"Saved Ewe-vs-I figure to: {iv_path}")
    plt.show()


if __name__ == "__main__":
    main()
