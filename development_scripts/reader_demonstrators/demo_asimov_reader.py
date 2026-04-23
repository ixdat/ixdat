import matplotlib.pyplot as plt

from ixdat import Measurement

DATASET_IDS = [
    "53ca3c5f-cd75-406c-8d99-a51c397988d7",
    "32d304fa-9439-4e97-801e-26cc59c67d37",
]


def _measurement_label(measurement):
    asimov_meta = (measurement.metadata or {}).get("asimov", {})
    return (
        asimov_meta.get("dataset_label")
        or asimov_meta.get("dataset_id")
        or measurement.name
    )


def main():
    measurements = [
        Measurement.read(dataset_id, reader="asimov") for dataset_id in DATASET_IDS
    ]
    fig, (ax_ec, ax_cv) = plt.subplots(2, 1, constrained_layout=True)
    ax_ec_right = ax_ec.twinx()

    line_styles = ("-", "--")
    u_color = "C0"
    j_color = "C3"
    for idx, meas in enumerate(measurements):
        label = _measurement_label(meas)
        linestyle = line_styles[idx % len(line_styles)]
        t_u, u = meas.grab(meas.U_name)
        ax_ec.plot(
            t_u,
            u,
            linestyle=linestyle,
            color=u_color,
            label=f"{label} U",
        )

        t_j, j = meas.grab(meas.J_name)
        ax_ec_right.plot(
            t_j,
            j,
            linestyle=linestyle,
            color=j_color,
            label=f"{label} J",
        )

    ax_ec.set_xlabel("time / [s]")
    if measurements:
        ax_ec.set_ylabel(measurements[0].U_name)
        ax_ec_right.set_ylabel(measurements[0].J_name)
    ax_ec.legend()
    ax_ec_right.legend()

    for idx, meas in enumerate(measurements):
        meas.as_cv().plot(
            ax=ax_cv,
            color=f"C{idx}",
            linestyle=line_styles[idx % len(line_styles)],
            label=_measurement_label(meas),
        )
    ax_cv.legend()

    plt.show()


if __name__ == "__main__":
    main()
