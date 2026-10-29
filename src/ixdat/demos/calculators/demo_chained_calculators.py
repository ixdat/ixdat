"""Demo of chaining calculators: MS quantification, background subtraction and EC
calibration applied together. Requires the demo data.

See https://github.com/ixdat/ixdat/issues/183

@author: Søren
"""

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.demos import get_demo_data_dir
from ixdat.exceptions import SeriesNotFoundError
from ixdat.techniques.ms import MSCalibration


def load_ecms(data_dir=None):
    """Return an ECMSMeasurement read from an ixdat .csv export."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(
        root / "ixdat_exports/trimarco2018_fig3_data.csv", reader="ixdat"
    )


EXAMPLES = {"ecms": load_ecms}


def main(show=True):
    results = {name: load() for name, load in EXAMPLES.items()}
    ecms = results["ecms"]

    try:
        ecms.grab("n_dot_H2", tspan=[100, 200])
    except SeriesNotFoundError as e:
        print("got error message, as intended:")
        print(e)

    # Chaining MS quantification and background subtraction
    ecms.add_calculator(MSCalibration("H2", "M2", F=0.21))
    t_H2, n_dot_H2 = ecms.grab("n_dot_H2", tspan=[100, 200])
    ecms.set_bg(mass_list=["M2"], tspan=[20, 30])
    n_dot_H2_bg_removed = ecms.grab_for_t("n_dot_H2", t=t_H2)
    n_dot_H2_raw = ecms.grab_for_t("n_dot_H2", t=t_H2, remove_background=False)

    M2_bg_removed = ecms.grab_for_t("M2", t=t_H2, remove_background=True)
    M2_raw = ecms.grab_for_t("M2", t=t_H2, remove_background=False)
    F_bg_implied = (M2_raw - M2_bg_removed) / (n_dot_H2_raw - n_dot_H2_bg_removed)
    assert F_bg_implied[0] == 0.21

    # "potential" is an alias for "raw_potential" when there are no calculators
    t, U_raw = ecms.grab("potential")
    print(f"max(U_raw) = {max(U_raw)}")

    ecms.calibrate(RE_vs_RHE=0.72)
    # ECMeasurement.calibrate() raises a warning and uses the last RE_vs_RHE:
    cal1a = ecms.calibrate(RE_vs_RHE=0.715)
    U = ecms.grab_for_t("potential", t=t)
    print(f"max(U) = {max(U)}")

    # ECMeasurement.calibrate() adds R_Ohm to the last calibration
    ecms.calibrate(R_Ohm=100)
    U_corr = ecms.grab_for_t("potential", t=t)
    print(f"max(U_corr) = {max(U_corr)}")

    # The "-raw" suffix ensures that no calculators are applied:
    U_raw_again = ecms.grab_for_t("potential-raw", t=t)
    print(f"max(U_raw_again) = {max(U_raw_again)}")
    U_again = ecms.grab_for_t("potential", calculator_list=[cal1a], t=t)
    print(f"max(U_again) = {max(U_again)}")

    assert max(U_again) == max(U)
    assert max(U_corr) != max(U)
    assert max(U_raw_again) == max(U_raw)
    assert max(U_raw) != max(U)

    if show:
        fig, ax = plt.subplots()
        ax.plot(t_H2, n_dot_H2, "b--", label="with background")
        ax.plot(t_H2, n_dot_H2_bg_removed, "b", label="background removed")
        ax.plot(t_H2, n_dot_H2_raw, "k--", label="raw")
        ax.set_xlabel("time / [s]")
        ax.set_ylabel("flux / [nmol/s]")
        ax.legend()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
