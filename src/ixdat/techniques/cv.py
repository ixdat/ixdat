import numpy as np
from .ec import ECMeasurement
from ..data_series import ValueSeries, TimeSeries
from ..exceptions import BuildError, SeriesNotFoundError
from .analysis_tools import (
    tspan_passing_through,
    calc_sharp_v_scan,
    find_signed_sections,
)
from ..plotters.ec_plotter import CVDiffPlotter
from ..plotters.plotting_tools import get_color_from_cmap, add_colorbar
from ..tools import deprecate


class CyclicVoltammogram(ECMeasurement):
    """Class for cyclic voltammetry measurements.

    Onto ECMeasurement, this adds:
    - a property `cycle` which is a ValueSeries on the same TimeSeries as potential,
    which counts cycles. "cycle" becomes the Measurement's `sel_str`. Indexing with
    integer or iterable selects according to `cycle`.
    - functions for quantitatively comparing cycles (like a stripping cycle, base cycle)
    - the default plot() is plot_vs_potential()
    """

    essential_series_names = ("t", "raw_potential", "raw_current", "cycle")
    selector_name = "cycle"

    series_constructors = ECMeasurement.series_constructors
    series_constructors["scan_rate"] = "_build_scan_rate"

    """Name of the default selector"""

    def __init__(self, *args, **kwargs):
        """Only reason to have an __init__ here is to set the default plot()"""
        super().__init__(*args, **kwargs)
        self.plot = self.plotter.plot_vs_potential  # gets the right docstrings! :D

        try:
            _ = self["cycle"]
        except SeriesNotFoundError:
            median_potential = 1 / 2 * (np.max(self.U) + np.min(self.U))
            self.redefine_cycle(start_potential=median_potential, redox=True)

        self.start_potential = None  # see `redefine_cycle`
        self.redox = None  # see `redefine_cycle`

    def __getitem__(self, key):
        """Given int list or slice key, return a CyclicVoltammogram with those cycles"""
        if isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            if step is None:
                step = 1
            key = list(range(start, stop, step))
        if isinstance(key, (int, list)):
            if isinstance(key, list) and not all([isinstance(i, int) for i in key]):
                print("can't get an item of type list unless all elements are int")
                print(f"you tried to get key = {key}.")
                raise AttributeError
            return self.select(key)
        return super().__getitem__(key)

    def redefine_cycle(self, start_potential=None, redox=None, N_points=5):
        """Build `cycle` which iterates when passing through start_potential

        Args:
            start_potential (float): The potential in [V] at which the cycle counter will
                iterate. If start_potential is not given, the cycle is just the
                `selector` inherited from ECMeasurement shifted to start at 0.
            redox (bool): True (or 1) for anodic, False (or 0) for cathodic. The
                direction in which the potential is scanning through start_potential to
                trigger an iteration of `cycle`.
            N_points (int): The number of consecutive points for which the potential
                needs to be above (redox=True) or below (redox=False) the
                start_potential for the new cycle to register.
        """
        self.start_potential = start_potential
        self.redox = redox
        if start_potential is None:
            old_cycle_series = self["cycle_number"]
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
            v = self.U
            if not redox:
                # easiest way to reverse directions is to use the same > < operators
                # but negate the arguments
                start_potential = -start_potential
                v = -v
            while n < N:
                # mask on remaining potential, True wherever behind the start potential:
                mask_behind = v[n:] < start_potential
                if True not in mask_behind:
                    # if the potenential doesn't go behind start potential again, then
                    # there are no more cycles
                    break
                else:
                    # the potential has to get behind the start potential for at least
                    # N_points data points before a new cycle can start.
                    n += np.argmax(mask_behind) + N_points

                # a mask on remaining potential, True wherever ahead of start potential:
                mask_in_front = v[n:] > start_potential
                if True not in mask_in_front:  # again, no more cycles.
                    break
                else:
                    # We've already been behind for N_points, so as soon as the
                    # potential gets ahead of the start_potential, a new cycle begins!
                    n += np.argmax(mask_in_front)
                c += 1
                cycle_vec[n:] = c  # and subsequent points increase in cycle number
                n += N_points  # have to be above start_potential for N_points
                # datapoints before getting behind it for this to count as a cycle.
            new_cycle_series = ValueSeries(
                name="cycle",
                unit_name="",
                data=cycle_vec,
                tseries=self.potential.tseries,
            )
        self.replace_series("cycle", new_cycle_series)

    def select_sweep(self, vspan, t_i=None):
        """Return the cut of the CV for which the potential is sweeping through vspan

        Args:
            vspan (iter of float): The range of self.potential for which to select data.
                Vspan defines the direction of the sweep. If vspan[0] < vspan[-1], an
                oxidative sweep is returned, i.e. one where potential is increasing.
                If vspan[-1] < vspan[0], a reductive sweep is returned.
            t_i (float): Optional. Time before which the sweep can't start
        """
        tspan = tspan_passing_through(
            t=self.t,
            v=self.U,
            vspan=vspan,
            t_i=t_i,
        )
        return self.cut(tspan=tspan)

    def integrate(self, item, tspan=None, vspan=None, ax=None):
        """Return the time integral of item while time in tspan or potential in vspan

        Args:
            item (str): The name of the ValueSeries to integrate
            tspan (iter of float): A time interval over which to integrate it
            vspan (iter of float): A potential interval over which to integrate it
        """
        if vspan:
            return self.select_sweep(
                vspan=vspan, t_i=tspan[0] if tspan else None
            ).integrate(item, ax=ax)
        return super().integrate(item, tspan, ax=ax)

    def _build_scan_rate(self, res_points=10):
        """The scan rate as a ValueSeries"""
        t, v = self.grab("potential")
        scan_rate_vec = calc_sharp_v_scan(t, v, res_points=res_points)
        scan_rate_series = ValueSeries(
            name="scan rate",
            unit_name="V/s",  # TODO: unit = potential.unit / potential.tseries.unit
            data=scan_rate_vec,
            tseries=self.potential.tseries,
        )
        return scan_rate_series

    @property
    @deprecate("0.1", "Use a look-up, i.e. `ec_meas['scan_rate']`, instead.", "0.3")
    def scan_rate(self):
        return self["scan_rate"]

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
            self["scan_rate"].data, x_res=v_scan_res, res_points=res_points
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
            vspan (iter of floats): The potential range in [V] to use for capacitance
        """
        sweep_1 = self.select_sweep(vspan)
        v_scan_1 = np.mean(sweep_1.grab("scan_rate")[1])  # [V/s]
        I_1 = np.mean(sweep_1.grab("raw_current")[1]) * 1e-3  # [mA] -> [A]

        sweep_2 = self.select_sweep([vspan[-1], vspan[0]], t_i=max(sweep_1.t + 1))
        v_scan_2 = np.mean(sweep_2.grab("scan_rate")[1])  # [V/s]
        I_2 = np.mean(sweep_2.grab("raw_current")[1]) * 1e-3  # [mA] -> [A]

        cap = 1 / 2 * (I_1 / v_scan_1 + I_2 / v_scan_2)  # [A] / [V/s] = [C/V] = [F]
        return cap

    def diff_with(self, other, v_list=None, cls=None, v_scan_res=0.001, res_points=10):
        """Return a CyclicVotammagramDiff of this CyclicVotammagram with another one

        Each anodic and cathodic sweep in other is lined up with a corresponding sweep
        in self. Each variable given in v_list (defaults to just "current") is
        interpolated onto self's potential and subtracted from self.

        Args:
            other (CyclicVoltammogram): The cyclic voltammogram to subtract from self.
            v_list (list of str): The names of the series to calculate a difference
                between self and other for (defaults to just "current").
            cls (ECMeasurement subclass): The class to return an object of. Defaults to
                CyclicVoltammogramDiff.
            v_scan_res (float): see :meth:`get_timed_sweeps`
            res_points (int):  see :meth:`get_timed_sweeps`
        """

        if not type(self) is CyclicVoltammogram:
            raise NotImplementedError(
                "CyclicVoltammogram.diff_with() is not implemented for "
                f"cyclic voltammograms of type {type(self)}"
            )

        vseries = self.potential
        tseries = vseries.tseries
        series_list = [tseries, self["raw_potential"], self["cycle"]]

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
                "Can only make diff of CyclicVoltammograms with same number of sweeps."
                f"{self!r} has {my_sweep_specs} and {other!r} has {others_sweep_specs}."
            )

        diff_values = {name: np.array([]) for name in v_list}
        t_diff = np.array([])

        for my_spec, other_spec in zip(my_sweep_specs, others_sweep_specs):
            sweep_type = my_spec[1]
            if not other_spec[1] == sweep_type:
                raise BuildError(
                    "Corresponding sweeps must be of same type when making diff."
                    f"Can't align {self!r}'s {my_spec} with {other!r}'s {other_spec}."
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

        cls = cls or CyclicVoltammogramDiff
        diff = cls.from_dict(diff_as_dict)
        # TODO: pass cv_compare_1 and cv_compare_2 to CyclicVoltammogramDiff as dicts
        diff.cv_compare_1 = self
        diff.cv_compare_2 = other
        return diff

    def plot_cycles(self, ax=None, cmap_name="jet"):
        """Plot the cycles on a color scale.

        Args:
            ax (mpl.Axis): The axes to plot on. A new one is made by default
            cmap_name (str): The name of the colormap to use. Defaults to "jet", which
                ranges from blue to red
        """
        cycle_numbers = set(self["cycle"].data)
        c_max = max(cycle_numbers)
        for c in cycle_numbers:
            color = get_color_from_cmap(c / c_max, cmap_name=cmap_name)
            ax = self[int(c)].plot(ax=ax, color=color)
        add_colorbar(
            ax, cmap_name, vmin=min(cycle_numbers), vmax=c_max, label="cycle number"
        )
        return ax


class CyclicVoltammagram(CyclicVoltammogram):

    # FIXME: decorating the class itself doesn't work because the callable returned
    #   by the decorator does not have the class methods. But this works fine.
    @deprecate("0.1", "Use `CyclicVoltammogram` instead ('o' replaces 'a').", "0.3")
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CyclicVoltammogramDiff(CyclicVoltammogram):

    default_plotter = CVDiffPlotter
    cv_compare_1 = None
    cv_compare_2 = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plot = self.plotter.plot
        self.plot_diff = self.plotter.plot_diff
        self.plotter = CVDiffPlotter(measurement=self)
