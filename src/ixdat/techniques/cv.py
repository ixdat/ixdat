from .ec import ECMeasurement


class CyclicVoltammagram(ECMeasurement):
    def plot(self, *args, **kwargs):
        """Default plot for cv is plot_vs_potential"""
        return self.plotter.plot_vs_potential(*args, **kwargs)
