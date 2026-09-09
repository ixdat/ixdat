"""Demonstrate the `tstamp_source` option of the OceanView reader.

OceanView writes a `Date:` line into the file header when it opens the file, and
an absolute timestamp on every spectral data row. The header is written a beat
before the first spectrum is actually acquired, so the two disagree by the time
it takes to open the file and start capturing.

`tstamp_source="header"` (the default) anchors the time axis to the header.
`tstamp_source="spectrum"` anchors it to the first spectrum's own timestamp,
which matters when aligning the optical data against another instrument.

Not every file can use "spectrum": older exports stamp their rows with the
1970-01-01 placeholder rather than a real date. The reader detects that and
refuses, rather than placing the measurement in 1970.

This script uses the two OceanView fixtures in test_data/oceanview_sec/ to show
both behaviours.
"""

from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt

from ixdat import Spectrum


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DATA_DIR = REPO_ROOT / "test_data" / "oceanview_sec"

# Rows carry real absolute timestamps, so both tstamp_source values work.
DATED_FILE = DATA_DIR / "spectrum_t_demo_QEP100221__0__18-51-43-086.txt"
# Rows carry the 1970-01-01 placeholder, so only "header" works.
PLACEHOLDER_FILE = DATA_DIR / "mini_oceanview__0__15-02-35-123.txt"


def show(tstamp):
    """Format a unix timestamp as an unambiguous UTC string."""
    return (
        datetime.fromtimestamp(tstamp, timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[
            :-3
        ]
        + " UTC"
    )


def compare_tstamp_sources():
    """Read one file both ways and report the difference."""
    from_header = Spectrum.read(DATED_FILE, reader="oceanview")
    from_spectrum = Spectrum.read(
        DATED_FILE, reader="oceanview", tstamp_source="spectrum"
    )

    shift = from_spectrum.tstamp - from_header.tstamp

    print(f"reading {DATED_FILE.name}")
    print(f"  tstamp_source='header'   -> tstamp = {show(from_header.tstamp)}")
    print(f"  tstamp_source='spectrum' -> tstamp = {show(from_spectrum.tstamp)}")
    print(f"  difference               : {shift:.6f} s")
    print()
    print("  this file's header states BST, and BST is one hour ahead of UTC,")
    print("  so a header of 18:51 BST is shown above as 17:51 UTC")
    print()
    print("  the header is written before acquisition starts, so 'spectrum'")
    print("  places the data later, where the spectrometer says it belongs")
    return from_header, from_spectrum


def show_tstamps_agree(spectrum_series, label):
    """ixdat holds the start time on the series and on the field's time axis.

    Both must carry the same value, or relative times get shifted when the data
    is accessed later.
    """
    series_tstamp = spectrum_series.tstamp
    field_tstamp = spectrum_series.field.axes_series[0].tstamp
    t_first = spectrum_series.field.axes_series[0].data[0]
    print(f"  {label}")
    print(f"    SpectrumSeries.tstamp : {show(series_tstamp)}")
    print(f"    Field time axis tstamp: {show(field_tstamp)}")
    print(f"    agree                 : {series_tstamp == field_tstamp}")
    print(f"    first relative time   : {t_first} s")


def demonstrate_placeholder_rejection():
    """Show that a 1970-stamped file is readable, but not with "spectrum"."""
    print(f"reading {PLACEHOLDER_FILE.name}")
    from_header = Spectrum.read(PLACEHOLDER_FILE, reader="oceanview")
    print(f"  tstamp_source='header'   -> tstamp = {show(from_header.tstamp)}")
    print("    the 1970 rows are still used for the relative time axis, where")
    print("    the placeholder date cancels out of the row-to-row differences")

    try:
        Spectrum.read(PLACEHOLDER_FILE, reader="oceanview", tstamp_source="spectrum")
    except ValueError as error:
        print("  tstamp_source='spectrum' -> ValueError:")
        for line in str(error).split(". "):
            print(f"    {line.strip()}")
    else:
        raise AssertionError("expected the placeholder date to be rejected")


def plot_the_offset(from_header, from_spectrum):
    """Plot both results on the same time axis, measured from the header
    timestamp."""
    wavelengths = from_header.field.axes_series[1].data
    index = len(wavelengths) // 2
    wavelength = wavelengths[index]

    fig, ax = plt.subplots()
    for series, label, style in (
        (from_header, "tstamp_source='header'", "o--"),
        (from_spectrum, "tstamp_source='spectrum'", "s-"),
    ):
        t_absolute = series.tstamp + series.field.axes_series[0].data
        ax.plot(
            t_absolute - from_header.tstamp,
            series.field.data[:, index],
            style,
            label=label,
        )
    ax.set_xlabel("time / [s] measured from the header timestamp")
    ax.set_ylabel(f"intensity at {wavelength:.0f} nm / [a.u.]")
    ax.legend()
    ax.set_title("the same spectra, anchored two ways")
    return ax


def main():
    print("=" * 70)
    print("comparing the two tstamp sources")
    print("=" * 70)
    from_header, from_spectrum = compare_tstamp_sources()

    print()
    print("=" * 70)
    print("both stored timestamps stay consistent")
    print("=" * 70)
    show_tstamps_agree(from_header, "tstamp_source='header'")
    show_tstamps_agree(from_spectrum, "tstamp_source='spectrum'")

    print()
    print("=" * 70)
    print("a file that cannot use the spectrum timestamp")
    print("=" * 70)
    demonstrate_placeholder_rejection()

    print()
    print("plotting the offset")
    plot_the_offset(from_header, from_spectrum)
    plt.show()


if __name__ == "__main__":
    main()
