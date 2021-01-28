"""This module contains loose functions and stuff useful for ixdat plotting."""


def color_axis(ax, color, lr="right", xy="y"):
    """Color the spine, ticks, and labels of an axis.

    ax (matplotlib.pyplot.axis): the axis to color (a part of)
    color (str): The color to color the axis.
    lr (str): whether to color the "left" spine or the "right". Defaults to "right".
    xy (str): whether to color the "x" axis or the "y". Defaults to "y".
    """
    ax.spines[lr].set_color(color)
    ax.tick_params(axis=xy, color=color)
    ax.tick_params(axis=xy, labelcolor=color)
    if xy == "y":
        ax.yaxis.label.set_color(color)
    if xy == "x":
        ax.xaxis.label.set_color(color)
