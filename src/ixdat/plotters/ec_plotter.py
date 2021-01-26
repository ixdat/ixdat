from matplotlib import pyplot as plt


class ECPlotter:
    def __init__(self, measurement=None):
        self.measurement = measurement

    def plot(self, measurement=None, tspan=None, ax=None):
        measurement = measurement or self.measurement
        t_v, v = measurement.get_potential(tspan=tspan)
        t_j, j = measurement.get_current(tspan=tspan)
        V_str = measurement.V_str
        J_str = measurement.J_str
        if not ax:
            fig, ax = plt.subplots()

        ax2 = ax.twinx()
        ax.plot(t_v, v, "k-", label=V_str)
        ax2.plot(t_j, j, "r-", label=J_str)
        ax.set_xlabel("time / [s]")
        ax.set_ylabel(V_str)
        ax2.set_ylabel(J_str)
