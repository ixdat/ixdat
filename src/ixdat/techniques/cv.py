import numpy as np
from .ec import ECMeasurement
from ..data_series import ValueSeries, TimeSeries
from ..exceptions import SeriesNotFoundError, BuildError
from .analysis_tools import (
    tspan_passing_through,
    calc_sharp_v_scan,
    find_signed_sections,
)


class CyclicVoltammagram(ECMeasurement):
    """Class for cyclic voltammatry measurements.

    Onto ECMeasurement, this adds:
    - a property `cycle` which is a ValueSeries on the same TimeSeries as potential,
    which counts cycles. "cycle" becomes the Measurement's `sel_str`. Indexing with
    integer or iterable selects according to `cycle`.
    - functions for quantitatively comparing cycles (like a stripping cycle, base cycle)
    - the default plot() is plot_vs_potential()
    """

    def __init__(self, *args, **kwargs):
        """Only reason to have an __init__ here is to set the default plot()"""
        super().__init__(*args, **kwargs)
        self.plot = self.plotter.plot_vs_potential  # gets the right docstrings! :D

    start_potential = None  # see `redefine_cycle`
    redox = None  # see `redefine_cycle`

    def __getitem__(self, key):
        """Given int list or slice key, return a CyclicVoltammagram with those cycles"""
        if type(key) is slice:
            start, stop, step = key.start, key.stop, key.step
            if step is None:
                step = 1
            key = list(range(start, stop, step))
        if type(key) in [int, list]:
            if type(key) is list and not all([type(i) is int for i in key]):
                print("can't get an item of type list unless all elements are int")
                print(f"you tried to get key = {key}.")
                raise AttributeError
            return self.select(key)
        try:
            return super().__getitem__(item=key)
        except SeriesNotFoundError:
            if key == "cycle":
                return self.cycle

    @property
    def cycle(self):
        """ValueSeries: the cycle number. The default selector. see `redefine_cycle`"""
        try:
            return self.selector
        except TypeError:
            # FIXME: This is what happens now when a single-cycle CyclicVoltammagram is
            #   saved and loaded.
            return ValueSeries(
                name="cycle",
                unit_name="",
                data=np.ones(self.t.shape),
                tseries=self.potential.tseries,
            )

    def redefine_cycle(self, start_potential=None, redox=None):
        """Build `cycle` which iterates when passing through start_potential

        Args:
            start_potential (float): The potential in [V] at which the cycle counter will
                iterate. If start_potential is not given, the cycle is just the
                `selector` inherited from ECMeasurement shifted to start at 0.
            redox (bool): True (or 1) for anodic, False (or 0) for cathodic. The
                direction in which the potential is scanning through start_potential to
                trigger an iteration of `cycle`.
        """
        self.start_potential = start_potential
        self.redox = redox
        if start_potential is None:
            old_cycle_series = self.cycle
            new_cycle_series = ValueSeries(
                name="cycle",
                unit_name=old_cycle_series.unit_name,
                data=old_cycle_series.data - min(old_cycle_series.data),
                tseries=old_cycle_series.tseries,
            )
        else:
            cycle_vec = np.zeros(self.t.shape)
            c = 0
            n = 0
            N = len(self.t)
            v = self.v
            if not redox:
                # easiest way to reverse directions is to use the same > < operators
                # but negate the arguments
                start_potential = -start_potential
                v = -v
            while n < N:
                mask_behind = v[n:] < start_potential
                if True not in mask_behind:
                    break
                else:
                    n += (
                        np.argmax(mask_behind) + 5
                    )  # have to be below V for 5 datapoints
                # print('point number on way up: ' + str(n)) # debugging

                mask_in_front = v[n:] > start_potential
                if True not in mask_in_front:
                    break
                else:
                    n += np.argmax(mask_in_front)
                c += 1
                cycle_vec[n:] = c  # and subsequent points increase in cycle number
                n += +5  # have to be above V for 5 datapoints
                # print('point number on way down: ' + str(n)) # debugging
            new_cycle_series = ValueSeries(
                name="cycle",
                unit_name="",
                data=cycle_vec,
                tseries=self.potential.tseries,
            )
        self["cycle"] = new_cycle_series
        self.sel_str = "cycle"
        return self.cycle

    def select_sweep(self, vspan, t_i=None):
        """Return a CyclicVoltammagram for while the potential is sweeping through vspan

        Args:
            vspan (iter of float): The range of self.potential for which to select data.
                Vspan defines the direction of the sweep. If vspan[0] < vspan[-1], an
                oxidative sweep is returned, i.e. one where potential is increasing.
                If vspan[-1] < vspan[0], a reductive sweep is returned.
            t_i (float): Optional. Time before which the sweep can't start.
        """
        tspan = tspan_passing_through(t=self.t, v=self.v, vspan=vspan, t_i=t_i,)
        return self.cut(tspan=tspan)

    def integrate(self, item, tspan=None, vspan=None, ax=None):
        """Return the time integral of item while time in tspan or potential in vspan

        item (str): The name of the ValueSeries to integrate
        tspan (iter of float): A time interval over which to integrate it
        vspan (iter of float): A potential interval over which to integrate it.
        """
        if vspan:
            return self.select_sweep(
                vspan=vspan, t_i=tspan[0] if tspan else None
            ).integrate(item, ax=ax)
        return super().integrate(item, tspan, ax=ax)

    @property
    def scan_rate(self, res_points=10):
        """The scan rate as a ValueSeries"""
        t, v = self.grab("potential")
        scan_rate_vec = calc_sharp_v_scan(t, v, res_points=res_points)
        scan_rate_series = ValueSeries(
            name="scan rate",
            unit_name="V/s",  # TODO: unit = potential.unit / potential.tseries.unit
            data=scan_rate_vec,
            tseries=self.potential.tseries,
        )
        # TODO: cache'ing, index accessibility
        return scan_rate_series

    def get_timed_sweeps(self, v_scan_res=5e-4, res_points=10):
        """Return list of [(tspan, type)] for all the potential sweeps in self.

        There are three types: "anodic" (positive scan rate), "cathodic" (negative scan
        rate), and "hold" (zero scan rate)

        Args:
            v_scan_res (float): The minimum scan rate considered significantly different
                than zero, in [V/s]. Defaults to 5e-4 V/s (0.5 mV/s). May need be higher
                for noisy potential, and lower for very low scan rates.
            res_points (int): The minimum number of points to be considered a sweep.
                During a sweep, a potential difference of at least `v_res` should be
                scanned through every `res_points` points.
        """
        t = self.t
        ec_sweep_types = {
            "positive": "anodic",
            "negative": "cathodic",
            "zero": "hold",
        }
        indexed_sweeps = find_signed_sections(
            self.scan_rate.data, x_res=v_scan_res, res_points=res_points
        )
        timed_sweeps = []
        for (i_start, i_finish), general_sweep_type in indexed_sweeps:
            timed_sweeps.append(
                ((t[i_start], t[i_finish]), ec_sweep_types[general_sweep_type])
            )
        return timed_sweeps

    def calc_capacitance(self, vspan):
        """Return the capacitance in [F], calculated by the first sweeps through vspan

        Args:
            vspan (iterable): The potential range in [V] to use for capacitance
        """
        sweep_1 = self.select_sweep(vspan)
        v_scan_1 = np.mean(sweep_1.grab("scan_rate")[1])  # [V/s]
        I_1 = np.mean(sweep_1.grab("raw_current")[1])  # [mA] -> [A]

        sweep_2 = self.select_sweep([vspan[-1], vspan[0]])
        v_scan_2 = np.mean(sweep_2.grab("scan_rate")[1])  # [V/s]
        I_2 = np.mean(sweep_2.grab("raw_current")[1]) * 1e-3  # [mA] -> [A]

        cap = 1/2 * (I_1 / v_scan_1 + I_2 / v_scan_2)  # [A] / [V/s] = [C/V] = [F]
        return cap

    def diff_with(self, other, v_list=None, cls=None, v_scan_res=0.001, res_points=10):
        """Return a CyclicVotammagramDiff of this CyclicVotammagram with another one

        Each anodic and cathodic sweep in other is lined up with a corresponding sweep
        in self. Each variable given in v_list (defaults to just "current") is
        interpolated onto self's potential and subtracted from self.

        Args:
            other (CyclicVoltammagram): The cyclic voltammagram to subtract from self.
            v_list (list of str): The names of the series to calculate a difference
                between self and other for (defaults to just "current").
            cls (ECMeasurement subclass): The class to return an object of. Defaults to
                CyclicVoltammagramDiff.
            v_scan_res (float): see CyclicVoltammagram.get_timed_sweeps()
            res_points (int):  see CyclicVoltammagram.get_timed_sweeps()
        """

        vseries = self.potential
        tseries = vseries.tseries
        series_list = [tseries, self.raw_potential, self.cycle]

        v_list = v_list or ["current", "raw_current"]
        if "potential" in v_list:
            raise BuildError(
                f"v_list={v_list} is invalid. 'potential' is used to interpolate."
            )

        my_sweep_specs = [
            spec
            for spec in self.get_timed_sweeps(
                v_scan_res=v_scan_res, res_points=res_points
            )
            if spec[1] in ["anodic", "cathodic"]
        ]
        others_sweep_specs = [
            spec
            for spec in other.get_timed_sweeps(
                v_scan_res=v_scan_res, res_points=res_points
            )
            if spec[1] in ["anodic", "cathodic"]
        ]
        if not len(my_sweep_specs) == len(others_sweep_specs):
            raise BuildError(
                "Can only make diff of CyclicVoltammagrams with same number of sweeps."
                f"{self} has {my_sweep_specs} and {other} has {others_sweep_specs}."
            )

        diff_values = {name: np.array([]) for name in v_list}
        t_diff = np.array([])

        for my_spec, other_spec in zip(my_sweep_specs, others_sweep_specs):
            sweep_type = my_spec[1]
            if not other_spec[1] == sweep_type:
                raise BuildError(
                    "Corresponding sweeps must be of same type when making diff."
                    f"Can't align {self}'s {my_spec} with {other}'s {other_spec}."
                )
            my_tspan = my_spec[0]
            other_tspan = other_spec[0]
            my_t, my_potential = self.grab(
                "potential", my_tspan, include_endpoints=False
            )
            t_diff = np.append(t_diff, my_t)
            other_t, other_potential = other.grab(
                "potential", other_tspan, include_endpoints=False
            )
            if sweep_type == "anodic":
                other_t_interp = np.interp(
                    np.sort(my_potential), np.sort(other_potential), other_t
                )
            elif sweep_type == "cathodic":
                other_t_interp = np.interp(
                    np.sort(-my_potential), np.sort(-other_potential), other_t
                )
            else:
                continue
            for name in v_list:
                my_v = self.grab_for_t(name, my_t)
                other_v = other.grab_for_t(name, other_t_interp)
                diff_v = my_v - other_v
                diff_values[name] = np.append(diff_values[name], diff_v)

        t_diff_series = TimeSeries(
            name="time/[s] for diffs", unit_name="s", data=t_diff, tstamp=self.tstamp
        )  # I think this is the same as self.potential.tseries

        series_list.append(t_diff_series)
        for name, data in diff_values.items():
            series_list.append(
                ValueSeries(
                    name=name,
                    unit_name=self[name].unit_name,
                    data=data,
                    tseries=t_diff_series,
                )
            )

        diff_as_dict = self.as_dict()
        del diff_as_dict["s_ids"]

        diff_as_dict["series_list"] = series_list
        diff_as_dict["raw_current_names"] = ("raw_current",)

        cls = cls or CyclicVoltammagramDiff
        diff = cls.from_dict(diff_as_dict)
        diff.cv_1 = self
        diff.cv_2 = other
        return diff


class CyclicVoltammagramDiff(CyclicVoltammagram):

    cv_1 = None
    cv_2 = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plot = self.plotter.plot
        self.plot_diff = self.plotter.plot_diff

    @property
    def plotter(self):
        """The default plotter for CyclicVoltammagramDiff is CVDiffPlotter"""
        if not self._plotter:
            from ..plotters.ec_plotter import CVDiffPlotter

            self._plotter = CVDiffPlotter(measurement=self)

        return self._plotter
