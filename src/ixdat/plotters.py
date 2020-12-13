from matplotlib import pyplot as plt


class ValuePlotter:
    def __init__(self, measurement=None):
        self.measurement = measurement

    def plot(self, *args, **kwargs):
        return self.plot_measurement(measurement=self.measurement, *args, **kwargs)

    def plot_measurement(self, measurement, v_list=None, tspan=None):
        fig, ax = plt.subplots()
        v_list = v_list or measurement.value_names

        for v_name in v_list:
            v, t = measurement.get_t_and_v(v_name, tspan=tspan)
            ax.plot(v, t, label=v_name)

        return ax
