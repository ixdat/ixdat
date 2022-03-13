"""Base class for plotters using matplotlib"""

from matplotlib import pyplot as plt
from matplotlib import gridspec


class MPLPlotter:
    """Base class for plotters based on matplotlib. Has methods for making mpl axes."""

    def new_ax(self, xlabel=None, ylabel=None):
        """Return a new matplotlib axis optionally with the given x and y labels"""
        fig, ax = plt.subplots()
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        return ax

    def new_two_panel_axes(self, n_bottom=1, n_top=1, emphasis="top"):
        """Return the axes handles for a bottom and top panel.

        TODO: maybe fix order of axes returned.
            see https://github.com/ixdat/ixdat/pull/30/files#r811198719

        Args:
            n_top (int): 1 for a single y-axis, 2 for left and right y-axes on top panel
            n_bottom (int): 1 for a single y-axis, 2 for left and right y-axes on bottom
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels

        Returns list of axes: top left, bottom left(, bottom right)(, top right)
        """
        self.new_ax()  # necessary to avoid deleting an open figure, I don't know why
        if emphasis == "top":
            gs = gridspec.GridSpec(5, 1)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:3, 0])]
            axes += [plt.subplot(gs[3:5, 0])]
        elif emphasis == "bottom":
            gs = gridspec.GridSpec(5, 1)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:2, 0])]
            axes += [plt.subplot(gs[2:5, 0])]
        else:
            gs = gridspec.GridSpec(6, 1)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:3, 0])]
            axes += [plt.subplot(gs[3:6, 0])]

        axes[0].xaxis.set_label_position("top")
        axes[0].tick_params(
            axis="x", top=True, bottom=False, labeltop=True, labelbottom=False
        )

        if n_bottom == 2:
            axes += [axes[1].twinx()]
        if n_top == 2:
            axes += [axes[0].twinx()]

        return axes
