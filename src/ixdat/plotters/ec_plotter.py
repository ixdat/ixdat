from matplotlib import pyplot as plt


class ECPlotter:
    def __init__(self, measurement=None):
        self.measurement = measurement

    def plot(self, measurement=None, tspan=None, V_str=None, J_str=None, axes=None):
        measurement = measurement or self.measurement
        V_str = V_str or measurement.V_str
        J_str = J_str or measurement.J_str
        t_v, v = measurement.get_t_and_v(V_str, tspan=tspan)
        t_j, j = measurement.get_t_and_v(J_str, tspan=tspan)
        if axes:
            ax1, ax2 = axes
        else:
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
        ax1.plot(t_v, v, "k-", label=V_str)
        ax2.plot(t_j, j, "r-", label=J_str)
        ax1.set_xlabel("time / [s]")
        ax1.set_ylabel(V_str)
        ax2.set_ylabel(J_str)
        return axes
