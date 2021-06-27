import numpy as np
from .ec import ECMeasurement
from ..data_series import ValueSeries
from ..exceptions import SeriesNotFoundError


class CyclicVoltammagram(ECMeasurement):
    """Class for cyclic voltammatry measurements.

    Onto ECMeasurement, this adds:
    - a property `cycle` which is a ValueSeries on the same TimeSeries as potential,
        which counts cycles. "cycle" becomes the Measurement's `sel_str`. Indexing with
        integer or iterable selects according to `cycle`.
    - functions for quantitatively comparing cycles (like a stripping cycle, base cycle)
    - the default plot() is plot_vs_potential()
    """

    sel_str = "cycle"
    """Name of the default selector"""

    def __init__(self, *args, **kwargs):
        """Only reason to have an __init__ here is to set the default plot()"""
        super().__init__(*args, **kwargs)
        self.plot = self.plotter.plot_vs_potential  # gets the right docstrings! :D

        self.start_potential = None  # see `redefine_cycle`
        self.redox = None  # see `redefine_cycle`

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
        return super().__getitem__(key)

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
        if "cycle" in self._cached_series:
            del [self._cached_series]
        self["cycle"] = new_cycle_series
