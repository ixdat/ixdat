"""Demo of the mass transport impulse response and deconvolution of EC-MS data.

Requires the demo data and the ``spectro_inlets_quantification`` (siq) package.

@author: awiniwar
"""

import matplotlib.pyplot as plt
import numpy as np

from ixdat import Measurement, plugins
from ixdat.calculators.ecms_calculators import MSCalResult, ECMSImpulseResponse
from ixdat.demos import get_demo_data_dir

# Start times of the O2 pulses. Accuracy to the second seems to be necessary.
PULSE_STARTS = [18172, 18472, 18772, 19072]
# Duration of the measured and modelled pulses, used to normalize the pulse area to 1
PULSE_LENGTH = 100
# End of the pulse, used for the baseline correction and the labels. It can be
# longer than PULSE_LENGTH if there are signal spikes which affect normalization.
PULSE_INF = 100


def load_o2_pulses(data_dir=None):
    """Return an MSMeasurement with a series of O2 pulses in a regular cell."""
    root = get_demo_data_dir(data_dir)
    return Measurement.read(root / "deconvolution/O2_pulses.csv", reader="ixdat")


def load_ca_dark(data_dir=None):
    """Return an ECMSMeasurement of a chronoamperometry on hematite in the dark."""
    root = get_demo_data_dir(data_dir)
    meas = Measurement.read(
        root / "deconvolution/to_deconvolute_data.csv", reader="ixdat"
    )
    meas.calibrate(A_el=0.35**2 * np.pi)
    return meas


EXAMPLES = {"o2_pulses": load_o2_pulses, "ca_dark": load_ca_dark}


def make_o2_calibration():
    """Return a siq Calculator with an O2 calibration, ready to use on measurements."""
    F_o2 = MSCalResult(
        name="O2", mol="O2", mass="M32", cal_type="ECMS calibration", F=0.11755
    ).to_siq()
    calibration = plugins.siq.Calculator(cal_list=[F_o2])
    calibration.set_quantifier(mol_list=["O2"], mass_list=["M32"], carrier="He")
    return calibration


def main(show=True):
    plugins.activate_siq()
    results = {name: load() for name, load in EXAMPLES.items()}
    pulses = results["o2_pulses"]
    ca_dark = results["ca_dark"]

    # Adding the calibration makes the flux of O2 available as "n_dot_O2"
    calibration = make_o2_calibration()
    pulses.add_calculator(calibration)

    # ---- Impulse response fitting of measured data, compared with a model
    if show:
        pulses.plot_measurement()
        fig, ax = plt.subplots()
        ax.set_xlabel("time / [s]")
        ax.set_ylabel("norm. intensity")

    for t_impulse in PULSE_STARTS:
        pulse = pulses.cut(tspan=[t_impulse - 10, t_impulse + PULSE_INF])
        # Shift the timestamp so that all pulses start at 0. The number is the start
        # of the full dataset, which was cut to make the file smaller.
        pulse.tstamp += 0.9355739 + t_impulse
        pulse_int = (
            pulse.integrate_flux(mol="O2", tspan=[0, PULSE_INF], tspan_bg=[-5, 0]) * 1e9
        )
        imp_resp = ECMSImpulseResponse.from_measurement(
            mol="O2",
            measurement=pulse,
            tspan=[-5, PULSE_LENGTH],
            tspan_bg=[[-5, 0], [PULSE_INF - 5, PULSE_INF]],
        )
        if show:
            ax.plot(
                imp_resp.t_kernel,
                imp_resp.kernel,
                marker="o",
                ms=3,
                markerfacecolor="w",
                linestyle="",
                label="n$_{O2}$" + " = {:.2f} nmol".format(pulse_int),
            )

    # Model impulse response. Vary the working distance (wd) by hand to fit the data.
    # With HOR in the same cell, a working distance of 230 um was measured.
    wd = 220e-6
    imp_resp_model = ECMSImpulseResponse.from_parameters(
        mol="O2",
        working_distance=wd,
        A_el=0.197,
        D=2.1e-9,  # optional if using siq
        H_v_cc=33,  # optional if using siq
        n_dot=None,  # optional if using siq
        carrier_gas=None,  # defaults to He
        gas_volume=1e-10,  # optional if using siq
        duration=PULSE_LENGTH,  # the same as the measured pulses, for normalization
        dt=0.1,
    )
    if show:
        ax.plot(
            imp_resp_model.t_kernel,
            imp_resp_model.kernel,
            linestyle="--",
            label="model = {:.0f} $\\mu$m".format(wd * 1e6),
        )
        ax.legend()

    # ---- Deconvolution of EC-MS data
    if show:
        ca_dark.plot(mass_list=["M32"], logplot=False)

    # A different working distance, as it is a different day and cell
    imp_resp_model = ECMSImpulseResponse.from_parameters(
        mol="O2",
        working_distance=180e-6,
        A_el=0.35**2 * np.pi,
        D=2.1e-9,  # optional if using siq
        H_v_cc=33,  # optional if using siq
        carrier_gas=None,  # defaults to He
        gas_volume=1e-10,  # optional if using siq
    )
    ca_dark.add_calculator(calibration)
    t_deconvoluted, v_deconvoluted = imp_resp_model.calc_deconvoluted_signal(
        mol="O2",
        measurement=ca_dark,
        tspan=[-10, 210],
        tspan_bg=[-10, -1],
    )
    if show:
        fig3, ax3 = plt.subplots()
        ax3.plot(t_deconvoluted, v_deconvoluted)

        # Deconvolute several tspans. Each tspan needs to include 0 if no t_zero is
        # given. This saves the plots under `name`.
        imp_resp_model.deconvolute_for_tspans(
            tspan_list=[[-10, 210], [208, 580]],
            t_zero_list=[6.6, 227],
            measurement=ca_dark,
            mol="O2",
            name="CA_hematite_dark_day1_deconvoluted_180um",
            export_data=False,
        )

    # Added as a calculator, the deconvolution makes the deconvoluted flux of its
    # molecule available as f"n_dot_{mol}-deconvoluted"
    ca_dark.add_calculator(imp_resp_model)
    t, n_dot_O2_decon = ca_dark.grab("n_dot_O2-deconvoluted", tspan_bg=[0, 10])

    if show:
        axes = ca_dark.plot(mol_list=["O2"], tspan_bg=[0, 10], logplot=False)
        axes[0].plot(t, n_dot_O2_decon, color="0.5", label="O2 deconvoluted")
        axes[0].legend()
        plt.show()
    return results


if __name__ == "__main__":
    results = main()
