import argparse

import matplotlib.pyplot as plt

from ixdat import Measurement


def demo(args):
    measurement = Measurement.read(
        args.dataset_id,
        reader="asimov",
        version=args.version,
        version_id=args.version_id,
        force_login=args.force_login,
    )
    cv = measurement.as_cv()

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, (ax_ec, ax_cv) = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
    fig.suptitle(measurement.name)

    measurement.plot(axes=[ax_ec, ax_ec.twinx()], U_color="#1f77b4", J_color="#d62728")
    cv.plot(ax=ax_cv, linewidth=1.5, color="#0f766e")

    if args.savefig:
        fig.savefig(args.savefig, dpi=220, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
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
        "--force-login",
        action="store_true",
        help="Force a fresh Keycloak device login",
    )
    parser.add_argument(
        "--savefig",
        default=None,
        help="Optional path to save figure image",
    )
    args = parser.parse_args()
    demo(args)
