import numpy as np
from .ec import ECMeasurement
from ..data_series import ValueSeries
from ..exceptions import SeriesNotFoundError


class CyclicVoltammagram(ECMeasurement):
    def plot(self, *args, **kwargs):
        """Default plot for cv is plot_vs_potential"""
        return self.plotter.plot_vs_potential(*args, **kwargs)

    def __getitem__(self, key):
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
        return self.selector

    def redefine_cycle(self, start_potential=None, redox=None):
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
            if redox in [0, -1, "red", "reduction"]:
                # easiest way to reverse directions is to use the same > < operators
                # but negate the arguments
                start_potential = -start_potential
                v = -v
            while n < N:
                mask_behind = v[n:] < start_potential
                if not True in mask_behind:
                    break
                else:
                    n += (
                        np.argmax(mask_behind) + 5
                    )  # have to be below V for 5 datapoints
                # print('point number on way up: ' + str(n)) # debugging

                mask_in_front = v[n:] > start_potential
                if not True in mask_in_front:
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
