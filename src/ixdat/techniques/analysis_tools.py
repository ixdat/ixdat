import numpy as np


def tspan_passing_through(t, v, vspan, direction=None, t_i=None, edge=None):

    t_i = t_i if t_i is not None else t[0] - 1
    edge = edge if edge is not None else np.abs(vspan[-1] - vspan[0]) / 100

    # define some things to generalize between anodic and cathodic
    if direction is None:
        direction = vspan[0] < vspan[-1]

    edge = edge if direction else -edge

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
                t > t_i, before(v, vspan[0] - edge)
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
