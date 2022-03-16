"""Miscellaneous tools for data analysis used in measurement techniques"""

import numpy as np
from scipy.optimize import minimize


def tspan_passing_through(t, v, vspan, direction=None, t_i=None, v_res=None):
    """Return the tspan corresponding to t when v first passes through vspan

    Args:
        t (np.array): independent varible data (usually time)
        v (np.array): dependent variable data
        vspan (iter of float): The range of v that we are interested in.
        direction (bool): Whether v should be increasing (True) or decreasing
            (False) as it passes through vspan. By default, the direction is
            defined by whether vspan is increasing or decreasing.
        t_i (float): The lowest value of t acceptable for tspan. Optional.
        v_res (float): The uncertainty or resolution of the v data. v must be
            at in or out of vspan by at least v_res to be considered in or out.
    """

    t_i = t_i if t_i is not None else t[0] - 1

    # define some things to generalize between anodic and cathodic
    if direction is None:
        direction = vspan[0] < vspan[-1]

    v_res = v_res if v_res is not None else np.abs(vspan[-1] - vspan[0]) / 100
    v_res = np.abs(v_res) if direction else -np.abs(v_res)

    def before(a, b):
        if direction:
            # before means more cathodic if we want the anodic sweep
            return a < b
        else:
            # and more anodic if we want the cathodic sweep
            return a > b

    if direction:
        # we start with the lower limit of V_span if we want the anodic sweep
        vspan = np.sort(np.array(vspan))
    else:
        # and with the upper limit of V_span if we want the cathodic sweep
        vspan = -np.sort(-np.array(vspan))

    t_before = t[
        np.argmax(
            np.logical_and(
                t > t_i, before(v, vspan[0] - v_res)
            )  # True if after t_i and comfortably out on start side
        )  # first index for which V is comfortably out on start side
    ]  # corresponding time
    t_just_before = t[
        np.argmax(
            np.logical_and(
                t > t_before, np.logical_not(before(v, vspan[0]))
            )  # True if after t_i and in on start side
        )
        - 1  # last index for which V is out on start side
    ]  # corresponding time
    i_start = np.argmax(np.logical_and(t > t_just_before, before(vspan[0], v)))
    # ^ first index of full sweep through range
    t_start = t[i_start]
    # ^ corresponding time
    i_finish = np.argmax(np.logical_and(t > t_start, before(vspan[1], v))) - 1
    # ^ last index of full sweep through range
    t_finish = t[i_finish]
    # ^ corresponding time
    return [t_start, t_finish]


def calc_sharp_v_scan(t, v, res_points=10):
    """Calculate the discontinuous rate of change of v with respect to t

    Args:
        t (np.array): the data of the independent variable, typically time
        v (np.array): the data of the dependent variable
        res_points (int): the resolution in data points, i.e. the spacing used in
            the slope equation v_scan = (v2 - v1) / (t2 - t1)
    """
    # the scan rate is dV/dt. This is a numerical calculation of dV/dt:
    v_behind = np.append(np.tile(v[0], res_points), v[:-res_points])
    v_ahead = np.append(v[res_points:], np.tile(v[-1], res_points))

    t_behind = np.append(np.tile(t[0], res_points), t[:-res_points])
    t_ahead = np.append(t[res_points:], np.tile(t[-1], res_points))

    v_scan_middle = (v_ahead - v_behind) / (t_ahead - t_behind)
    # ^ this is "softened" at the anodic and cathodic turns.

    # We can "sharpen" it by selectively looking ahead and behind:
    v_scan_behind = (v - v_behind) / (t - t_behind)
    v_scan_ahead = (v_ahead - v) / (t_ahead - t)

    # but this gives problems right at the beginning, so set those to zeros
    v_scan_behind[:res_points] = np.zeros(res_points)
    v_scan_ahead[-res_points:] = np.zeros(res_points)

    # now sharpen the scan rate!
    v_scan = v_scan_middle
    mask_use_ahead = np.logical_and(
        np.abs(v_scan_ahead) > np.abs(v_scan),
        np.abs(v_scan_ahead) > np.abs(v_scan_behind),
    )
    v_scan[mask_use_ahead] = v_scan_ahead[mask_use_ahead]

    mask_use_behind = np.logical_and(
        np.abs(v_scan_behind) > np.abs(v_scan),
        np.abs(v_scan_behind) > np.abs(v_scan_ahead),
    )
    v_scan[mask_use_behind] = v_scan_behind[mask_use_behind]

    return v_scan


def find_signed_sections(x, x_res=0.001, res_points=10):
    """Return list of tuples ((i_start, i_finish), section_type) describing the vector x

    `i_start` and `i_finish` are indexes in x defining where sections start and end.
    `section_type` can be "positive" (x>0), "negative" (x<0) or "zero" (x~0).

    Args:
        x (np array): The data as a vector
        x_res (float): The minimum value in x to be considered different from zero,
            i.e. the uncertainty or resolution of the data
        res_points (int): The minimum number of consecutive data points in x that must
            have the same sign (or be ~0) to constitute a section of the data.
    """
    mask_negative = x < -x_res
    mask_positive = x > x_res
    mask_zero = abs(x) < x_res

    section_types = ["negative", "positive", "zero"]
    the_masks = [mask_negative, mask_positive, mask_zero]

    for mask in the_masks:
        mask[-2] = False
        mask[-1] = True
    N = len(x)
    i_start = 0
    i_finish = 0
    n_sweep = 0

    the_next_starts = [np.argmax(mask) for mask in the_masks]
    section_id = int(np.argmin(the_next_starts))

    sections = []
    while i_start < N - 1:
        I_out = np.argmin(the_masks[section_id][i_finish:])
        the_next_start = i_finish + I_out + res_points

        try:
            I_in_again = np.argmax(the_masks[section_id][the_next_start:])
        except ValueError:
            the_next_starts[section_id] = N
        else:
            the_next_starts[section_id] = the_next_start + I_in_again
            # ^ and add it.

        next_section_id = int(np.argmin(the_next_starts))
        i_finish = the_next_starts[next_section_id]

        if next_section_id != section_id:
            sections.append(((i_start, i_finish), section_types[section_id]))
            section_id = next_section_id
            n_sweep += 1
            i_start = i_finish
        else:
            i_start += res_points

    return sections


def calc_t_using_scan_rate(v, dvdt):
    """Return a numpy array describing the time corresponding to v given scan rate dvdt

    This is useful for data sets where time is missing. It depends on another value
    having a constant absolute rate of change (such as electrode potential in cyclic
    voltammatry).
    It uses the `calc_sharp_v_scan` algorithm to match the scan rate implied by the
    timevector returned with the given scan rate.
    Args:
        v (np array): The value
        dvdt (float): The scan rate in units of v's unit per second
    Returns:
        np array: t, the time vector corresponding to v
    """

    def error(t_tot):
        t = np.linspace(0, t_tot[0], v.size)
        dvdt_calc = np.abs(calc_sharp_v_scan(t, v))
        error = np.sum(dvdt_calc**2 - dvdt**2)
        return error

    t_total_guess = (max(v) - min(v)) / dvdt
    result = minimize(error, np.array(t_total_guess))

    t_total = result.x[0]
    return np.linspace(0, t_total, v.size)
