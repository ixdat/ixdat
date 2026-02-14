import argparse
import matplotlib.pyplot as plt

from ixdat import Measurement
from ixdat.exceptions import SeriesNotFoundError


def _measurement_label(measurement):
    asimov_meta = (measurement.metadata or {}).get("asimov", {})
    return (
        asimov_meta.get("dataset_label")
        or asimov_meta.get("dataset_id")
        or measurement.name
    )


def demo(args):
    dataset_ids = (
        args.dataset_id if isinstance(args.dataset_id, list) else [args.dataset_id]
    )
    if len(dataset_ids) == 1:
        measurement = Measurement.read(
            dataset_ids[0],
            reader="asimov",
            version=args.version,
            version_id=args.version_id,
            force_login=args.force_login,
        )
        measurements_to_plot = [measurement]
    else:
        measurement = Measurement.read_set(
            file_list=dataset_ids,
            reader="asimov",
            version=args.version,
            version_id=args.version_id,
            force_login=args.force_login,
        )
        measurements_to_plot = measurement.component_measurements or [measurement]

    fig, (ax_ec, ax_cv) = plt.subplots(2, 1, constrained_layout=True)
    ax_ec_right = ax_ec.twinx()

    line_styles = ("-", "--")
    u_color = "C0"
    j_color = "C3"
    for idx, meas in enumerate(measurements_to_plot):
        label = _measurement_label(meas)
        linestyle = line_styles[idx % len(line_styles)]
        try:
            t_u, u = meas.grab(meas.U_name)
        except SeriesNotFoundError:
            pass
        else:
            ax_ec.plot(
                t_u,
                u,
                linestyle=linestyle,
                color=u_color,
                label=f"{label} U",
            )

        try:
            t_j, j = meas.grab(meas.J_name)
        except SeriesNotFoundError:
            pass
        else:
            ax_ec_right.plot(
                t_j,
                j,
                linestyle=linestyle,
                color=j_color,
                label=f"{label} J",
            )

    ax_ec.set_xlabel("time / [s]")
    if measurements_to_plot:
        ax_ec.set_ylabel(measurements_to_plot[0].U_name)
        ax_ec_right.set_ylabel(measurements_to_plot[0].J_name)
    ax_ec.legend()
    ax_ec_right.legend()

    for idx, meas in enumerate(measurements_to_plot):
        meas.as_cv().plot(
            ax=ax_cv,
            color=f"C{idx}",
            linestyle=line_styles[idx % len(line_styles)],
            label=_measurement_label(meas),
        )
    ax_cv.legend()

    if args.savefig:
        fig.savefig(args.savefig, dpi=220, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read and plot one or more Asimov datasets using ixdat."
    )
    parser.add_argument(
        "--dataset-id",
        nargs="+",
        required=True,
        help="One or more Asimov dataset UUIDs",
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
