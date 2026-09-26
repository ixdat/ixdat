import numpy as np
from .ec import ECMeasurement
from ..data_series import ValueSeries, TimeSeries
from ..exceptions import BuildError, SeriesNotFoundError
from ..calculators.scan_rate_tools import (
    tspan_passing_through,
    find_signed_sections,
    calc_sharp_v_scan
)
from ..plotters import CVDiffPlotter, get_color_from_cmap, add_colorbar
from ..tools import deprecate

import warnings



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

    built_in_calculator_types = ECMeasurement.built_in_calculator_types + [
        "scan_rate_calculator"
    ]

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

    def redefine_cycle(
        self, start_potential=None, redox=None, N_points=5, turning_point=False, N_sep=10, res_points=None, rate_threshold=0.0005, N_points_threshold=1
    ):
        """Build `cycle` which iterates when passing through start_potential, or a turning point

        Args:
            start_potential (float): The potential in [V] at which the cycle counter will
                iterate. If start_potential is not given, the cycle is just the
                `selector` inherited from ECMeasurement shifted to start at 0.
            redox (bool): True (or 1) for anodic, False (or 0) for cathodic. If selecting via start_potential only, 
                this is the direction in which the potential is scanning through start_potential to
                trigger an iteration of `cycle`. For turning_points, a new cycle will iterate if the scan rate after the turning point
                matches the redox direction. If redox is neither true nor false and turning_points is 
                true, the cycle selector will iterate at every turning point. This can then defacto be used
                to select anodic or cathodic sweeps. If turning_points is false, a redox direction must be selected.
            N_points (int): If turning_point is False, this is the number of consecutive 
                points for which the potential needs to be above (redox=True) or below 
                (redox=False) the start_potential for the new cycle to register. If turning_point
                is True, this is the number of consecutive points after a detected turning
                point for which the scan rate turning point must be greater than the rate_threshold 
                (redox=True) or lower than the negative of the rate_threshold (redox=False), 
                or the same sign as the point directly after turning point
            turning_point (bool): If True, define cycles using changes in
                the direction of the potential sweep instead of a fixed
                potential.
            N_sep (int) : The number of indices by which the turning points must be separated. 10 by default.
            res_points (int): to be passed to calc_sharp_v. 10 by default in calc_sharp_v The resolution in data points,
                i.e. the spacing used in the slope equation v_scan = (v2 - v1) / (t2 - t1)
            rate_threshold (float) : the threshold magnitude in scan rate to accept for turning point
                validity. 0.5 mV/s by default.
            N_points_threshold (float) : Must be a value between 0 and 1. If data is noisy, this can be
                invoked to soften the requirement of N_points so that a fraction of them must meet the 
                requirements. For example, if 8 of 10 points meet the requirement but there are 2 noisy
                data points, the turning point will still be accepted.            
        """
        self.start_potential = start_potential
        self.redox = redox
        if turning_point:

            v = self.U
            N = len(v)
            time = self.t
            
            # define cycle vector
            cycle_vec = np.full(N, np.nan)
            
            #Find point where potential first crosses start potential, if given as an argument as well as turning_point
            if start_potential is not None:
                crossing = np.where(
                    (
                            (v[:-1] < start_potential) &
                            (v[1:] >= start_potential)
                    )
                    |
                    (
                            (v[:-1] > start_potential) &
                            (v[1:] <= start_potential)
                    )
                )[0]
                if len(crossing) == 0:
                    raise ValueError(
                        f"No crossing of "
                        f"{start_potential} V found."
                    )
                start_idx = crossing[0] + 1

                #define cycle 0
                cycle_vec[:start_idx] = 0
            else:
                start_idx = 0
                
            
            scan_rate=calc_sharp_v_scan(time, v, res_points=res_points)
            # Potential holds are likely to be noisy and oscillate around zero. Use rate_threshold to naively find areas of potential holds
            # as a first check.
            sign = np.zeros_like(scan_rate)
            sign[scan_rate>rate_threshold] = 1 
            sign[scan_rate<-rate_threshold] = -1
            
            if redox is True:
                # Negative -> positive
                turning_indices = np.where((sign[:-1] <= 0) & (sign[1:] > 0))[0]

            elif redox is False:
                # Positive -> negative
                turning_indices = np.where((sign[:-1] >= 0) & (sign[1:] < 0))[0]

            else:
                # Either direction
                warnings.warn("No redox direction selected. Cycles will iterate at every turning point.")
                turning_indices = np.where(sign[:-1] != sign[1:])[0]
                
            #Each cycle should have only 1 turning index. If there are multiple, this is due to noise.
            #This section checks that the following N_points all have either a positive (for anodic sweep)
            #or negative (for cathodic sweep) sign. If redox is not given, the code check that all the signs
            #for the indices directly after the turning point are the same.
            #The N_sep variable is used to check that the valid turning points are not closer together than
            #the user defined separation of N_sep
            if len(turning_indices) >= 1:
                valid_indices = []
                sign_checked_indices = []
                separation_clusters = []

                for idx in turning_indices:
                    signcheck=False
                    # Don't accept turning points before the chosen start
                    if idx < start_idx:
                        continue
                    window_end = idx +1 + N_points
                    next_points = scan_rate[idx+1:window_end]
                    if len(next_points) > 0:
                        if redox is True:
                            same_sign = next_points > rate_threshold
                        elif redox is False:
                            same_sign = next_points < -rate_threshold
                        else:
                            same_sign = np.sign(next_points) == sign[idx+1]
                        
                        if same_sign is not None:
                            if N_points_threshold<0 or N_points_threshold>1:
                                raise ValueError("N_points_threshold must be value between 0 and 1")
                            fraction_correct=np.mean(same_sign)
                            if fraction_correct>=N_points_threshold:
                                #valid_indices.append(idx)
                                signcheck=True
                        if signcheck:
                            sign_checked_indices.append(idx)
                if len(sign_checked_indices)==0:
                    raise ValueError("No valid turning points found")
                # Only loop over turning points that passed sign check to check for separation. Cluster indices that are too close together.
                # First point is included in the first cluster by default
                cluster = [sign_checked_indices[0]]
                for idx in sign_checked_indices[1:]:
                    if idx-cluster[-1]<N_sep:
                        cluster.append(idx)
                    else:
                        separation_clusters.append(cluster)
                        cluster = [idx]
                #Get the final cluster
                separation_clusters.append(cluster)
                
                #Choose the best turning point in this cluster. The best turning point is the one with the smallest scan rate.
                for cluster in separation_clusters:
                    valid_indices.append(cluster[np.argmin(scan_rate[cluster])])             
                
            turning_indices = np.asarray(valid_indices)


            for c, idx in enumerate(turning_indices, start=start_idx):
                cycle_vec[idx:] = c

            new_cycle_series = ValueSeries(
                name="cycle",
                unit_name="",
                data=cycle_vec,
                tseries=self.potential.tseries,
            )

        elif start_potential is None:
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
            if redox==None:
                raise ValueError("Redox direction must be selected if utilising potential as cycle indicator")
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

    @property
    @deprecate("0.1", "Use a look-up, i.e. `ec_meas['scan_rate']`, instead.", "0.3.1")
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

    def as_ec(self):
        """Convert self to an EC measurement"""
        from .ec import ECMeasurement

        ec_as_dict = self.as_dict()
        ec_as_dict["technique"] = "EC"
        # Note, this works perfectly! All needed information is in self_as_dict :)
        return ECMeasurement.from_dict(ec_as_dict)


class CyclicVoltammagram(CyclicVoltammogram):
    # FIXME: decorating the class itself doesn't work because the callable returned
    #   by the decorator does not have the class methods. But this works fine.
    @deprecate("0.1", "Use `CyclicVoltammogram` instead ('o' replaces 'a').", "0.3.1")
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
