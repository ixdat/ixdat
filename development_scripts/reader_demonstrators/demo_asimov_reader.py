import matplotlib.pyplot as plt

from ixdat import Measurement
from ixdat.exceptions import BuildError, SeriesNotFoundError, TechniqueError

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


def read_measurements(dataset_ids):
    return [
        Measurement.read(
            dataset_id,
            reader="asimov",
        )
        for dataset_id in dataset_ids
    ]


def resolve_common_series_key(measurements, candidate_keys, quantity_name):
    for key in candidate_keys:
        try:
            for measurement in measurements:
                measurement.grab(key)
        except SeriesNotFoundError:
            continue
        return key

    raise ValueError(
        f"These Asimov datasets cannot be plotted together with this demo. "
        f"No shared {quantity_name} lookup key was found. Tried: "
        + ", ".join(repr(key) for key in candidate_keys)
    )


def prepare_measurements_for_demo(measurements):
    if not measurements:
        raise ValueError("No Asimov measurements were loaded for plotting.")

    reference = measurements[0]
    reference_u_name = resolve_common_series_key(
        measurements,
        ["potential", reference.U_name, "raw_potential", reference.E_name],
        "voltage",
    )
    reference_j_name = resolve_common_series_key(
        measurements,
        ["current", reference.J_name, "raw_current", reference.I_name],
        "current",
    )
    cv_measurements = []
    incompatible = []

    for measurement in measurements:
        label = _measurement_label(measurement)
        issues = []

        try:
            t_u, u = measurement.grab(reference_u_name)
        except SeriesNotFoundError:
            issues.append(f"missing voltage series for key {reference_u_name!r}")
        else:
            if len(t_u) == 0 or len(u) == 0:
                issues.append(f"empty voltage series for key {reference_u_name!r}")

        try:
            t_j, j = measurement.grab(reference_j_name)
        except SeriesNotFoundError:
            issues.append(f"missing current series for key {reference_j_name!r}")
        else:
            if len(t_j) == 0 or len(j) == 0:
                issues.append(f"empty current series for key {reference_j_name!r}")

        try:
            cv_measurement = measurement.as_cv()
        except (BuildError, TechniqueError, SeriesNotFoundError, ValueError) as exc:
            issues.append(f"cannot convert to CV: {type(exc).__name__}: {exc}")
        else:
            cv_measurements.append(cv_measurement)

        if issues:
            incompatible.append(f"{label}: " + "; ".join(issues))

    if incompatible:
        raise ValueError(
            "These Asimov datasets cannot be plotted together with this demo. "
            "The overlay expects measurements that share the same voltage/current "
            "series names, have non-empty voltage/current data, and can all be "
            "converted with as_cv().\n- " + "\n- ".join(incompatible)
        )

    return reference_u_name, reference_j_name, cv_measurements


def plot_measurements(measurements_to_plot):
    reference_u_name, reference_j_name, cv_measurements = prepare_measurements_for_demo(
        measurements_to_plot
    )

    fig, (ax_ec, ax_cv) = plt.subplots(2, 1, constrained_layout=True)
    ax_ec_right = ax_ec.twinx()

    line_styles = ("-", "--")
    u_color = "C0"
    j_color = "C3"
    for idx, meas in enumerate(measurements_to_plot):
        label = _measurement_label(meas)
        linestyle = line_styles[idx % len(line_styles)]
        try:
            t_u, u = meas.grab(reference_u_name)
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
            t_j, j = meas.grab(reference_j_name)
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
        ax_ec.set_ylabel(reference_u_name)
        ax_ec_right.set_ylabel(reference_j_name)
    ax_ec.legend()
    ax_ec_right.legend()

    for idx, cv_measurement in enumerate(cv_measurements):
        cv_measurement.plot(
            ax=ax_cv,
            color=f"C{idx}",
            linestyle=line_styles[idx % len(line_styles)],
            label=_measurement_label(measurements_to_plot[idx]),
        )
    ax_cv.legend()

    plt.show()


def main():
    plot_measurements(read_measurements(DATASET_IDS))


if __name__ == "__main__":
    main()
