"""This module contains loose functions and stuff useful for ixdat plotting."""

import matplotlib as mpl
from matplotlib import pyplot as plt


def color_axis(ax, color, lr="right", xy="y"):
    """Color the spine, ticks, and labels of an axis.

    ax (matplotlib.pyplot.axis): the axis to color (a part of)
    color (str): The color to color the axis.
    lr (str): whether to color the "left" spine or the "right". Defaults to "right".
    xy (str): whether to color the "x" axis or the "y". Defaults to "y".
    """
    ax.spines[lr].set_color(color)
    ax.tick_params(axis=xy, color=color, labelcolor=color)
    if xy == "y":
        ax.yaxis.label.set_color(color)
    if xy == "x":
        ax.xaxis.label.set_color(color)


def get_color_from_cmap(x, cmap_name):
    """Return the color as a 4-tuple, given a value btwn 0 and 1, and color map name

    Args:
        x (float): Value between 0 and 1 defining a location on a matplotlib colormap
        cmap_name (str): The name of a matplotlib colormap. Popular ones include
            "inferno" (0=black --> 1=yellow) and "jet" (0=blue --> 1=red)
            See https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
    """

    cmap = mpl.cm.get_cmap(cmap_name)

    rgba = cmap(x)
    return rgba


def add_colorbar(ax, cmap_name, vmin, vmax, label="intensity"):
    """Add a colorbar to a matplotlib axis

    ax (matplotlib Axis): The axis (with data plotted on it) to add a colorbar for
    cmap_name (str): The name of the color map
    vmin (float): The minimum value in the plotted data
    vmax (float): The maximum value in the plotted data
    label (str): A label for the color bar (the name of the value represented by color)
    """

    cmap = mpl.cm.get_cmap(cmap_name)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = plt.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        use_gridspec=True,
        anchor=(0.75, 0),
    )
    cb.set_label(label)
