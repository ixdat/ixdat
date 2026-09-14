"""Tests for Matplotlib plotter layout helpers."""

from matplotlib import pyplot as plt

from ixdat.plotters.base_mpl_plotter import MPLPlotter


def assert_offset_text_is_inside(axis):
    """Assert that an axis displays its scientific scale inside its frame."""
    offset_text = axis.yaxis.get_offset_text()
    text_bounds = offset_text.get_window_extent(axis.figure.canvas.get_renderer())

    assert offset_text.get_text()
    assert text_bounds.y1 <= axis.bbox.y1
    assert text_bounds.y0 >= axis.bbox.y0


def test_three_panel_axes_share_x_without_a_placeholder_axis():
    """The three-panel helper creates only its returned axes."""
    axes = MPLPlotter().new_three_panel_axes(
        n_middle=2,
        n_bottom=2,
        interactive=False,
    )
    primary_axes = axes[:3]
    returned_axes = [axis for axis in axes if axis is not None]
    figure = primary_axes[0].figure

    try:
        assert figure.axes == returned_axes
        assert all(
            primary_axes[0].get_shared_x_axes().joined(primary_axes[0], axis)
            for axis in primary_axes[1:]
        )

        figure.canvas.draw()
        top_ticks = primary_axes[0].xaxis.get_major_ticks()
        middle_ticks = primary_axes[1].xaxis.get_major_ticks()
        bottom_ticks = primary_axes[2].xaxis.get_major_ticks()

        assert all(not tick.label1.get_visible() for tick in top_ticks)
        assert all(not tick.label2.get_visible() for tick in top_ticks)
        assert all(not tick.label1.get_visible() for tick in middle_ticks)
        assert all(not tick.label2.get_visible() for tick in middle_ticks)
        assert all(tick.label1.get_visible() for tick in bottom_ticks)
        assert all(not tick.label2.get_visible() for tick in bottom_ticks)

        primary_axes[0].set_xlabel("time / [s]")
        primary_axes[1].set_xlabel("time / [s]")
        assert not primary_axes[0].xaxis.label.get_visible()
        assert not primary_axes[1].xaxis.label.get_visible()

        for axis in (primary_axes[1], primary_axes[2], axes[4], axes[5]):
            axis.plot([0, 1], [1e-9, 2e-9])
        figure.canvas.draw()

        for axis in (primary_axes[1], primary_axes[2], axes[4], axes[5]):
            assert_offset_text_is_inside(axis)
    finally:
        plt.close(figure)


def test_two_panel_lower_axes_place_scale_text_inside():
    """Two panels share one bottom x-axis and keep lower scale text visible."""
    axes = MPLPlotter().new_two_panel_axes(n_bottom=2, interactive=False)
    figure = axes[0].figure

    try:
        assert axes[0].get_shared_x_axes().joined(axes[0], axes[1])
        axes[0].set_xlabel("time / [s]")
        figure.canvas.draw()

        top_ticks = axes[0].xaxis.get_major_ticks()
        bottom_ticks = axes[1].xaxis.get_major_ticks()
        assert all(not tick.label1.get_visible() for tick in top_ticks)
        assert all(not tick.label2.get_visible() for tick in top_ticks)
        assert all(tick.label1.get_visible() for tick in bottom_ticks)
        assert all(not tick.label2.get_visible() for tick in bottom_ticks)
        assert not axes[0].xaxis.label.get_visible()

        for axis in (axes[1], axes[3]):
            axis.plot([0, 1], [1e-9, 2e-9])
        figure.canvas.draw()

        for axis in (axes[1], axes[3]):
            assert_offset_text_is_inside(axis)
    finally:
        plt.close(figure)
