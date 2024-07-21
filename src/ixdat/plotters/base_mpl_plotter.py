"""Base class for plotters using matplotlib"""

from collections import defaultdict

from matplotlib import pyplot as plt
from matplotlib import gridspec


class MPLPlotter:
    """Base class for plotters based on matplotlib. Has methods for making mpl axes."""

    def __init__(self):
        # Instantiate data holders for dynamic range selection
        self._dynamically_added_objects = defaultdict(list)
        self._axis_for_range_selection = set()
        self._selected_range = {"left": None, "right": None}

    def new_ax(self, xlabel=None, ylabel=None, interactive=True):
        """Return a new matplotlib axis optionally with the given x and y labels

        Args:
            xlabel (str): The label to apply to the x-axis
            ylabel (str): The label to apply to the y-axis
            interactive (bool): Whether to activate interactive range selection (default
                True)

        """

        fig, ax = plt.subplots()
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)

        # Add the axis to those we perform range selection on and connect mouse events
        if interactive:
            self._axis_for_range_selection.add(ax)
            fig.canvas.mpl_connect("button_press_event", self.onclick)

        return ax

    def new_two_panel_axes(self, n_bottom=1, n_top=1, emphasis="top", interactive=True):
        """Return the axes handles for a bottom and top panel.

        Args:
            n_top (int): 1 for a single y-axis, 2 for left and right y-axes on top panel
            n_bottom (int): 1 for a single y-axis, 2 for left and right y-axes on bottom
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels
            interactive (bool): Whether to activate interactive range selection (default
                True)

        Returns list of axes: top left, bottom left(, top right, bottom right)
        """
        # Necessary to avoid deleting an open figure:
        fig = plt.figure()

        if emphasis == "top":
            gs = gridspec.GridSpec(5, 1, figure=fig)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:3, 0])]
            axes += [plt.subplot(gs[3:5, 0])]
        elif emphasis == "bottom":
            gs = gridspec.GridSpec(5, 1, figure=fig)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:2, 0])]
            axes += [plt.subplot(gs[2:5, 0])]
        else:
            gs = gridspec.GridSpec(6, 1, figure=fig)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:3, 0])]
            axes += [plt.subplot(gs[3:6, 0])]

        if interactive:
            self._axis_for_range_selection = set(axes)

        axes[0].xaxis.set_label_position("top")
        axes[0].tick_params(
            axis="x", top=True, bottom=False, labeltop=True, labelbottom=False
        )

        if n_bottom == 2 or n_top == 2:
            axes += [None, None]
        if n_top == 2:
            axes[2] = axes[0].twinx()
        if n_bottom == 2:
            axes[3] = axes[1].twinx()

        return axes

    def new_three_panel_axes(self, n_bottom=1, n_middle=1, n_top=1, interactive=True):
        """Return the axes handles for a bottom, middle, and top panel.

        Args:
            n_top (int): 1 for a single y-axis, 2 for left and right y-axes on top panel
            n_middle (int): 1 for a single y-axis, 2 for left and right y-axes on middle
            n_bottom (int): 1 for a single y-axis, 2 for left and right y-axes on bottom
            interactive (bool): Whether to activate interactive range selection (default
                True)

        Returns list of axes: top left, middle left, bottom left(, top right, middle
            right, bottom right)
        """
        # Necessary to avoid deleting an open figure, I don't know why
        self.new_ax(interactive=interactive)

        gs = gridspec.GridSpec(12, 1)
        # gs.update(hspace=0.025)
        axes = [plt.subplot(gs[0:4, 0])]
        axes += [plt.subplot(gs[4:8, 0])]
        axes += [plt.subplot(gs[8:12, 0])]

        if interactive:
            self._axis_for_range_selection = set(axes)

        axes[0].xaxis.set_label_position("top")
        axes[0].tick_params(
            axis="x", top=True, bottom=False, labeltop=True, labelbottom=False
        )
        axes[1].tick_params(
            axis="x", top=True, bottom=True, labeltop=False, labelbottom=False
        )

        if n_bottom == 2 or n_middle == 2 or n_top == 2:
            axes += [None, None, None]
        if n_top == 2:
            axes[3] = axes[0].twinx()
        if n_middle == 2:
            axes[4] = axes[1].twinx()
        if n_bottom == 2:
            axes[5] = axes[2].twinx()

        return axes

    def onclick(self, event):
        """Place range markers in plot"""
        # Don't place markers if outside the plotted area
        if event.xdata is None or event.ydata is None:
            return

        # Clear the previous marker line of this type (left/right)
        for line in self._dynamically_added_objects.pop(event.button, []):
            line.remove()

        # Just remove the marker on double-clicks
        if event.dblclick:
            self._selected_range[event.button.name.lower()] = None
            plt.draw()
            return

        # Add the new marker line
        for ax in self._axis_for_range_selection:
            ylim = ax.get_ylim()
            self._dynamically_added_objects[event.button] += ax.plot(
                [event.xdata] * 2,
                ylim,
                color="black",
                linewidth=0.2,
            )
            ax.set_ylim(ylim)

        # Add to recorded limits and print
        self._selected_range[event.button.name.lower()] = event.xdata
        if (
            self._selected_range["left"] is not None
            and self._selected_range["right"] is not None
        ):
            # When we have both left and right selection, extract the axis type and form
            # a nice range name
            extracted_xlabel = ""
            for ax in self._axis_for_range_selection:
                extracted_xlabel = ax.get_xlabel()
                if extracted_xlabel != "":
                    break

            if "time" in extracted_xlabel:
                range_name = "tspan"
            else:
                range_name = "xspan"

            # Print span and span size
            span_size = abs(self._selected_range["right"] - self._selected_range["left"])
            print(
                f"{range_name}={list(sorted(self._selected_range.values()))}"
                f"   span_size={span_size}"
            )
        else:
            # Print the one added selector
            for side, value in self._selected_range.items():
                if value is not None:
                    print(f"{side}={value}")

        plt.draw()
