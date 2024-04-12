"""This module contains loose functions and stuff useful for ixdat plotting."""

import numpy as np
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


def smooth_vector(y, n_points):
    """Return copy of the vector `y` smoothed by a running average of `n_points`"""
    # extend the vector on each side to avoid edge effects. The total extension is
    #   n_points long, split between start and finish.
    start_pad = int(np.floor(n_points))
    end_pad = n_points - start_pad
    y_extended = np.append(
        np.append(y[0] * np.ones((start_pad,)), y), y[-1] * np.ones((end_pad,))
    )
    # Note, the use of cumsum is faster than convolve, according to the answer here:
    #   https://stackoverflow.com/a/34387987
    cumsum_vec = np.cumsum(y_extended)  # put a 0 at the front of y
    y_smooth = (cumsum_vec[n_points:] - cumsum_vec[:-n_points]) / n_points
    return y_smooth


def calc_linear_background(t, y, tspans):
    """Return a copy of the vector `y` that interpolates linearly between tspans

    The vector `y - calc_linear_background(t, y, tspans)` will go to zero at the times
    on `t` specified by `tspans`

    Args:
        t (numpy Array): time
        y (numpy Array): the value to calculate a background to
        tspans (list of timespans): The times to interpolate the background between
    """
    t_bg_list = []
    y_bg_list = []
    for tspan in tspans:
        mask = np.logical_and(tspan[0] < t, t < tspan[-1])
        if True not in mask:
            continue
        t_bg_list.append(t[mask].mean())
        y_bg_list.append(y[mask].mean())
    return np.interp(t, t_bg_list, y_bg_list)


def get_indeces_and_times(
    t_vec,
    dt=None,
    t_list=None,
    dn=None,
    index_list=None,
):
    """Given a time vector, find the indeces and times of a desired subset

    Specify one of dt, t_list, dn, or index_list

    Args:
        t_vec (numpy array): The time vector
        dt (float): time interval between spectra to plot, [s]. The
            first spectrum and those taken at times closest to each
            integer multiple of dt after are plotted.
        t_list (list of float): List of times for which to plot the
            spectrum, [s]. The closest spectrum to each time in the
            list is plotted.
        dn (int): number of spectra between plotted spectra
        index_list (list of int): List of indeces of spectra to plot
    Returns Tuple:
        list of int: The indeces of the desired points
        list of float: the times of the desired points

    """
    if not dt and not t_list and not dn and not index_list:
        raise ValueError("must specify one of: dt, t_list, dn, index_list")
    if dt or t_list:
        if dt:
            num = int((t_vec[-1] - t_vec[0]) / dt)
            t_list = [t_vec[0] + n * dt for n in range(num)]
        index_list = []
        for t in t_list:
            index = int(np.argmin(np.abs(t_vec - t)))
            index_list.append(index)
    elif dn:
        index_list = list(range(0, len(t_vec), dn))
    if not t_list:
        t_list = [t_vec[i] for i in index_list]
    return index_list, t_list
