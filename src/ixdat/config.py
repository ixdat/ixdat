"""Configuration variables including standard data directory

This module also defines user-specific options, such as whether to use plugins

The object `plugins` created here gives direct access to these options.
Example useage:

```
import ixdat

ixdat.config.plugins.use_si_quant = True  # use the spectro_inlets_quantification package
```

See `help(ixdat.options.plugins)` for information.
"""
import datetime
import warnings
from pathlib import Path
import numpy as np
from .tools import deprecate


class _Config:
    """
    Attributes:
        standard_data_directory (Path): the directory in which to save by default with
            the default database backend, which is to save files in a directory.
            ixdat will make the directory if it does not exist.
        standard_metadata_suffix (str): The file ext. for JSON format metadata files
        standard_data_suffix (str): The file extension for numpy.save format data files
        timestamp_string_format (str): A format string for datetime.datetime.strftime.
            Defaults to ixdats custom datetime format: 22E18 14:34:55
        timezone (datetime.timezone): The timezone timestamps should use when formatted
            to string. Defaults to the current local timestamp.
    """

    def __init__(self):
        self.standard_metadata_suffix = ".ix"
        self.standard_data_suffix = ".ix.npy"
        self.standard_ixdat_directory = Path.home() / ".ixdat"
        self.standard_data_directory = self.standard_ixdat_directory / "projects"
        self.default_project_name = "test"
        self.timestamp_string_format = "native"
        self.timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo

    @property
    def ixdat_temp_dir(self):
        temp_dir = self.standard_data_directory / "temp"
        if not temp_dir.exists():
            temp_dir.mkdir(parents=True)
        return temp_dir


config = _Config()


def prompt_for_permission(prompt):
    yn = input(prompt + "\nEnter y for yes or anything else for no.")
    return yn in ["y", "yes", "Y", "Yes"]


class _PluginOptions:
    """A class for activating plugins and giving access to items from plugins

    This class has only one instance, initiated on import of ixdat. Use it by
    `from ixdat.config import plugins`

    These packages need to be separately installed.

    Packages
    --------

    spectro_inlets_quantification
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    - use_si_quant (bool): Read-only. If this is True, ixdat uses the
        `spectro_inlets_quantification` package. This changes the behaviour of some
        methods in `MSMeasurement` and inheriting classes. Defaults to False.
    - activate_si_quant(): Sets use_si_quant to True and initializes the quant. package
        classes and the siq `Calculator` are then available as `siq.<class_name>`
    - deactivate_si_quant(): Sets use_si_quant to False.
    - si_quant (_QuantDeps). This is where the needed imports from the external
        quantification package go. See the docstring of _QuantDeps

    cinfdata
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    - activate_cinfdata(): sets up a database as `cinfdata` The database has to be
        connected using `cinfdata.connect(...)`
    """

    def __init__(self):
        self._use_siq = False
        self.siq = _SIQ()
        self.cinfdata = _CinfData()

    @property
    def use_siq(self):
        return self._use_siq

    def activate_siq(self):
        """Changes mass spec methods to use spectro_inlets_quantification"""
        self._use_siq = True
        self.siq.populate()

    def deactivate_siq(self):
        """Changes mass spec methods to use ixdat's native MS quantification"""
        self._use_siq = False

    def activate_cinfdata(self):
        self.cinfdata.activate()

    @deprecate(
        last_supported_release="0.2.6",
        update_message="`siq` is the universal abreviation "
        "for `spectro_inlets_quantification`. "
        "Thus, `ixdat.plugins.activate_si_quant` "
        "is now `ixdat.plugins.activate_siq`",
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def activate_si_quant(self):
        return self.activate_siq()

    @property
    @deprecate(
        last_supported_release="0.2.6",
        update_message="`siq` is the universal abreviation "
        "for `spectro_inlets_quantification`. "
        "Thus, `ixdat.plugins.use_si_quant` "
        "is now `ixdat.plugins.use_siq`",
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def use_si_quant(self):
        return self.use_siq

    @property
    @deprecate(
        last_supported_release="0.2.6",
        update_message="`siq` is the universal abreviation "
        "for `spectro_inlets_quantification`. "
        "Thus, `ixdat.plugins.si_quant` is now `ixdat.plugins.siq`",
        hard_deprecation_release="0.3.0",
        remove_release="1.0.0",
    )
    def si_quant(self):
        return self.siq


class _SIQ:
    """Class storing items of `spectro_inlets_quantification`.

    This class has one instance, which is an attribute of `plugins`. To print this
    docstring, you would type:
    ```
    from ixdat.config import plugins

    plugins.use_si_quant = True  # Activates plugins.quant.
    help(plugins.si_quant)  #gives information on the si_quant package
    ```

    The attributes of this class are `None` until the property `plugins.use_si_quant` is
    set to True, triggering their population (activating quant).

    Once activated, the attributes of `plugins.quant` are:
    - `Chip`: Class describing the MS inlet. More powerful than ixdat's `MSInlet`
    - `Molecule`: Class with data about molecules relevant to (EC-)MS quantification
    - `CalPoint`: Class with data about an MS calibration experiment
    - `Calibration`: Class for storing, visualizing, and using multiple CalPoints.
    - `Quantifier`: Class for using a Calibration to quantify MS data
    - `quant_config`: The config object of the external quantification package

    `plugins.quant` also has
    - `QUANT_DIRECTORY`: A property for getting and setting the data directory used by
      the external quantification package.
    """

    def __init__(self):
        self.Chip = None
        self.Molecule = None
        self.Calibration = None
        self.CalPoint = None
        self.Quantifier = None
        self.quant_config = None
        self.Calculator = None
        self._QUANT_DIRECTORY = None

    def populate(self):
        from spectro_inlets_quantification.config import Config
        from spectro_inlets_quantification.medium import Medium
        from spectro_inlets_quantification.molecule import Molecule
        from spectro_inlets_quantification.chip import Chip
        from spectro_inlets_quantification.calibration import CalPoint, Calibration
        from spectro_inlets_quantification.quantifier import Quantifier
        from .measurements import Calculator
        from .data_series import TimeSeries, ValueSeries
        from .techniques.ms import MSCalculator
        from .techniques.ec_ms import ECMSCalculator
        from .exceptions import QuantificationError

        self.quant_config = Config()
        self.medium = Medium()
        self.Molecule = Molecule
        self.Chip = Chip
        self.CalPoint = CalPoint
        self.Calibration = Calibration
        self.Quantifier = Quantifier

        def native_to_siq_method(native_method):
            """A decorator easing the writing of the siqCalculator class"""

            def siq_method(*args, **kwargs):
                # `chip`s go by the more general name `inlet`s in native ixdat. (The
                # two are usually interchangable as the primary method is in both cases
                # called `calc_n_dot`.
                # Make sure to actually initiate and use a siq Chip if you want to
                # benifit from the more powerful capillary flux calculator in
                # the calibration methods.
                chip = kwargs.pop("chip", None)
                if chip:
                    kwargs["inlet"] = chip
                plugins.deactivate_siq()  # to suppress warning
                result = native_method(*args, **kwargs)
                plugins.activate_siq()
                return result.to_siq()

            return siq_method

        class siqCalculator(Calibration, Calculator):
            """A class to function as an ixdat calculator with siq functionality"""

            # FIXME: Should this somehow be in calclators/ms_calculators.py
            # FIXME: Should this be defined within the siq package?

            def __init__(self, *args, **kwargs):
                self.measurement = kwargs.pop("measurement", None)
                Calibration.__init__(self, *args, **kwargs)
                self.quantifier = None

            def set_quantifier(self, **kwargs):
                """Set the `spectro_inlets_quantification` quantifier.

                The Quantifier is an object with the method `calc_n_dot`, which takes a
                dictionary of signals or signal vectors in [A] and return a dictionary of
                molecular fluxes in [mol/s].
                The quantifier typically does this by solving the linear equations of
                S_M = sum_i ( F_M^i * n_dot^i )
                Where n_dot^i is the flux to the vacuum chamber of molecule i in [mol/s], S_M
                is the signal at mass M in [A], and F_M^i is the *sensitivity factor* of molecule
                i at mass M.
                The quantifier thus needs access to a set of sensitivity factors.

                The quantifier can be built in this method (avoiding explicit import of the
                `spectro_inlets_quantification` package) by providing the sensitivity factors
                in the form of a `siqCalculator` (which can be obtained from e.g.
                plugins.siq.Calculator.multicomp_gas_flux_cal) and the specification of which ones to
                use by `mol_list` and `mass_list`.
                The quantifier will always use all the masses in `mass_list` to solve for the
                flux of all the mols in `mol_list`.

                The argument `carrier` is required by some quantifiers but only used if
                partial pressures before the MS inlet are required (`quantifier.calc_pp`)

                Quantification is only as accurate as your sensitivity factors!

                Args:
                    mol_list (list of str): The list of molecules to use in flux calculations.
                       These should all be represented in the Calibration. If not provided,
                       we'll use all the mols in the Calibration.
                    mass_list (list of str): The list of masses to use in flux calculations.
                       These should all be represented in the Calibration. If not provided,
                       we'll use all the masses in the Calibration.
                    carrier (optional, str): The carrier gas in the experiment. Defaults to "He".

                See `help(plugins.siq.quantifier.__init__` for other keyword args
                """
                self.quantifier = Quantifier(calibration=self, **kwargs)

            @property
            def available_series_names(self):
                if not self.quantifier:
                    warnings.warn(
                        "a siqCalculator has no available series"
                        "until its quantifier is set using `set_quantifier`"
                    )
                    return set()
                return set([f"n_dot_{mol}" for mol in self.quantifier.mol_list])

            gas_flux_calibration = native_to_siq_method(
                MSCalculator.gas_flux_calibration
            )
            gas_flux_calibration_curve = native_to_siq_method(
                MSCalculator.gas_flux_calibration_curve
            )
            ecms_calibration = native_to_siq_method(ECMSCalculator.ecms_calibration)
            ecms_calibration_curve = native_to_siq_method(
                ECMSCalculator.ecms_calibration_curve
            )

            @classmethod
            def siq_multicomp_gas_flux_calibration(
                cls,
                measurement,
                mol_list,
                mass_list,
                gas,
                tspan,
                gas_bg=None,
                tspan_bg=None,
                chip=None,
            ):
                """Calibration of multiple components of a calibration gas simultaneously

                Uses a matrix equation and the reference spectra in the molecule data files.

                The results are only as accurate as the reference spectrum used. For this reason,
                this method is a last resort and it is recommended *not* to use a multicomponent
                calibration gas. Instead, get a separate calibration gas for each molecule to
                be calibrated.

                Here is an explanation of the math used in this method:

                The fundamental matrix equation is:
                  S_vec = F_mat @ n_dot_vec
                Elementwise, this is:
                 S_M = sum_i ( F^i_M * n_dot^i )
                Rewrite to show that sensitivity factors follow each molecule's spectrum:
                 S_M = sum_i (F_weight_i * spectrum^i_M * n_dot^i)
                And regroup the parts that only depend on the molecule (^i):
                 S_M = sum_i (spectrum^i_M * (F_weight^i * n_dot^i))
                 S_M = sum_i (spectrum^i_M * sensitivity_flux^i)
                Change back into a matrix equation, and solve it:
                 S_vec = spectrum_mat @ sensitivity_flux_vec
                 sensitivity_flux_vec = spectrum_mat^-1 @ S_vec   # eq. 1
                Ungroup the part we grouped before (the "sensitivity_flux"):
                 F_weight^i = sensitivity_flux^i / n_dot^i        # eq. 2
                And, in the end, each sensitivity factor is:
                 F_M^i = F_weight^i * spectrum^i_M                # eq. 3

                Equations 1, 2, and 3 are implemented in the code of this method.

                Args:
                    mol_list (list of str): List of the names of the molecules to calibrate
                    mass_list (list of str): List of the masses to calibrate
                    gas (Gas, dict, or str): Composition of the calibration gas, e.g.
                       {"Ar": 0.95, "H2": 0.05} for 5% H2 in Ar
                    tspan (Timespan): Timespan during which the calibration gas is in the chip
                    gas_bg (Gas, dict, or str): Composition of the background gas
                    tspan_bg (Timespan): Timespan during which the background gas is in the chip
                    chip (Chip, optional): object describing the MS capillary, if different than
                       the standard chip in the MS quantification package

                Returns Calibration: An object from `spectro_inlets_quantification`,
                   representing all the calibration results from the calibration.
                """
                if not plugins.use_siq:
                    raise QuantificationError(
                        "`MSMeasurement.siq_multicomp_gas_flux_calibration` "
                        "only works when using `spectro_inlets_quantification` "
                        "(`ixdat.plugins.activate_siq()`). "
                    )

                chip = chip or Chip()
                chip.gas = gas
                flux = chip.calc_n_dot()

                chip_bg = chip or Chip()
                chip_bg.gas = gas_bg
                flux_bg = chip_bg.calc_n_dot()

                delta_flux_list = []
                for mol in mol_list:
                    delta_flux = flux.get(mol, 0) - flux_bg.get(mol, 0)
                    delta_flux_list.append(delta_flux)
                delta_flux_vec = np.array(delta_flux_list)

                delta_signal_list = []
                for mass in mass_list:
                    S = self.grab_signal(mass, tspan=tspan)[1].mean()
                    if tspan_bg:
                        S_bg = self.grab_signal(mass, tspan=tspan_bg)[1].mean()
                    else:
                        S_bg = 0
                    delta_S = S - S_bg
                    delta_signal_list.append(delta_S)
                delta_signal_vec = np.array(delta_signal_list)

                spectrum_vec_list = []
                for mol in mol_list:
                    spectrum = chip.gas.mdict[mol].norm_spectrum
                    spectrum_vec = np.array(
                        [spectrum.get(mass, 0) for mass in mass_list]
                    )
                    spectrum_vec_list.append(spectrum_vec)
                spectrum_mat = np.stack(spectrum_vec_list).transpose()

                inverse_spectrum_mat = np.linalg.inv(spectrum_mat)
                sensitivity_flux_vec = inverse_spectrum_mat @ delta_signal_vec  # eq. 1
                F_weight_vec = sensitivity_flux_vec / delta_flux_vec  # eq. 2

                cal_list = []
                for i, mol in enumerate(mol_list):
                    for M, mass in enumerate(mass_list):
                        F = F_weight_vec[i] * spectrum_mat[M, i]  # eq. 3
                        if F:
                            cal = CalPoint(
                                mol=mol,
                                mass=mass,
                                F=F,
                                F_type="capillary",
                                date=self.yyMdd,
                            )
                            cal_list.append(cal)

                return cls(cal_list=cal_list, measurement=measurement)

            def grab_fluxes(
                self,
                measurement=None,
                tspan=None,
                tspan_bg=None,
                remove_background=None,
                include_endpoints=None,
            ):
                """Return a time vector and a dictionary with all the quantified fluxes

                Args:
                    tspan (list): Timespan for which the signal is returned.
                    tspan_bg (list): Timespan that corresponds to the background signal.
                        If not given, no background is subtracted.
                    remove_background (bool): Whether to remove a pre-set background if
                        available. Defaults to True..
                    include_endpoints (bool): Whether to interpolate for tspan[0] and
                        tspan[-1]
                """
                measurement = measurement or self.measurement
                if not self.quantifier:
                    raise QuantificationError(
                        "A siqCalculator object has to have its `Quantifier` defined"
                        "before it can calculate fluxes. Use method `set_quantifier`."
                    )
                sm = self.quantifier.sm
                signals = {}
                t = None
                for mass in sm.mass_list:
                    if t is None:
                        t, S = measurement.grab(
                            mass,
                            tspan=tspan,
                            tspan_bg=tspan_bg,
                            remove_background=remove_background,
                            include_endpoints=include_endpoints,
                        )
                    else:
                        S = measurement.grab_for_t(
                            mass,
                            t=t,
                            tspan_bg=tspan_bg,
                            remove_background=remove_background,
                        )
                    signals[mass] = S
                return sm.calc_n_dot(signals=signals)

            def calculate_series(self, key, measurement=None):
                measurement = measurement or self.measurement

                t, n_dots = self.grab_siq_fluxes()

                return ValueSeries(
                    name="n_dot_{key}",
                    data=n_dots[key],
                    unit="mol/s",
                    tseries=TimeSeries(
                        name="t", data=t, unit="s", tstamp=measurement.tstamp
                    ),
                )

        self.Calculator = siqCalculator

    @property
    def QUANT_DIRECTORY(self):
        if not self._QUANT_DIRECTORY:
            self._QUANT_DIRECTORY = (
                config.standard_ixdat_directory / "plugin_data/ms_quant"
            )
            if not self._QUANT_DIRECTORY.exists():
                self._QUANT_DIRECTORY.mkdir(parents=True)
        return self._QUANT_DIRECTORY

    @QUANT_DIRECTORY.setter
    def QUANT_DIRECTORY(self, quant_directory):
        self._QUANT_DIRECTORY = quant_directory
        self.quant_config.aux_data_directory = self.QUANT_DIRECTORY


class _CinfData:
    """Class implement direct database read access using external module Cinfdata"""

    def __init__(
        self,
    ):
        self.cinfdata = None
        self._managed_cinfdata_object = None
        self._context_manager_kwargs = None

    def activate(self):
        from cinfdata import Cinfdata

        self.DB = Cinfdata

    def connect(self, setup_name=None, grouping_column=None):
        """setup_name (str): The name of the table inside the database
        grouping_column (str): Either the 'timestamp'/'comment' or 'Comment' column"""
        return self.DB(setup_name=setup_name, grouping_column=grouping_column)

    def __call__(self, **kwargs):
        """**kwargs: setup_name (str) and grouping_column (str) (see connect())"""
        self._context_manager_kwargs = kwargs
        return self

    def __enter__(self):
        if self._managed_cinfdata_object:
            self._context_manager_kwargs = None
            raise RuntimeError(
                "Using the cinfdata plugin as a context manager more "
                "than once at the same time is not supported"
            )

        self._managed_cinfdata_object = self.connect(**self._context_manager_kwargs)
        return self._managed_cinfdata_object

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._managed_cinfdata_object.connection.close()
        self._managed_cinfdata_object = None
        self._context_manager_kwargs = None


plugins = _PluginOptions()
