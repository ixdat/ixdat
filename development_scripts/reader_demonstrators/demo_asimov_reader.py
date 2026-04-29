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
    for idx, meas in enumerate(measurements):
        linestyle = line_styles[idx % len(line_styles)]
        meas.plot(
            axes=[ax_ec, ax_ec_right],
            U_color="C0",
            J_color="C3",
            linestyle=linestyle,
        )

    for idx, meas in enumerate(measurements):
        meas.plot_vs_potential(
            ax=ax_cv,
            color=f"C{idx}",
            linestyle=line_styles[idx % len(line_styles)],
            label=_measurement_label(meas),
        )
    ax_cv.legend()

    plt.show()


if __name__ == "__main__":
    main()
