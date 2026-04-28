"""Demo script for the SPEC .dat reader (SpecDATReader)

Test data: FeFoil_Fe_exafs.dat — Fe foil at Fe K-edge (~7.1 keV),
3 repeat transmission scans from SOLEIL synchrotron.

Column assignment for this beamline:
  I0  = mon_4      (incident monitor, smooth across energy)
  It  = ion_1_2    (transmitted beam, shows Fe K-edge as absorption dip)
  µt  = -log(It/I0)

Note: raw count rates give It > I0 (different detector efficiencies), so µt
has a constant negative offset. The edge step and EXAFS oscillations are
physically correct — EXAFS analysis software normalises by the edge step.

Demonstrates:
  - MultiSpectrum (all columns, no technique)
  - Single-scan and multi-scan XAS
  - Scan selection by number
  - Multi-scan averaging
  - Metadata, tstamp, duration stored on the returned object
  - Helpful ValueError when I0/It columns are not specified
"""

from pathlib import Path
import datetime
import matplotlib.pyplot as plt
from ixdat import Spectrum

data_file = Path(__file__).parents[2] / "test_data/spec/FeFoil_Fe_exafs.dat"

I0_col = "mon_4"
It_col = "ion_1_2"

# =============================================================================
# 1. Inspect what the reader stores on the returned object
# =============================================================================

xas = Spectrum.read(
    data_file,
    reader="spec",
    technique="XAS",
    y_name=It_col,
    ref_name=I0_col,
)

print("=== SpecDATReader — stored attributes ===")
print(f"  name      : {xas.name}")
print(f"  technique : {xas.technique}")
print(f"  tstamp    : {xas.tstamp}  " f"({datetime.datetime.fromtimestamp(xas.tstamp)})")
print(f"  duration  : {xas.duration} s  " f"({xas.duration / 60:.1f} min)")
print(f"  reader    : {type(xas.reader).__name__}")
print()
meta = xas.metadata
print(f"  metadata['scan_number']  : {meta['scan_number']}")
print(f"  metadata['scan_command'] : {meta['scan_command']}")
print(f"  metadata['count_time_ms']: {meta['count_time_ms']} ms")
print(f"  metadata['beamline_file']: {meta['beamline_file']}")
print(f"  metadata['comment']      : {meta['comment']}")
print("  metadata['umi']          :")
for line in meta["umi"]:
    print(f"    {line}")
print(f"  metadata['motor_positions'] ({len(meta['motor_positions'])} motors):")
for k, v in list(meta["motor_positions"].items())[:6]:
    print(f"    {k}: {v}")
print("    ...")
print()

# =============================================================================
# 2. Demonstrate ValueError when required columns are omitted
# =============================================================================

print("=== ValueError when y_name / ref_name not given ===")
try:
    Spectrum.read(data_file, reader="spec", technique="XAS")
except ValueError as e:
    print(f"  {e}")
print()

# =============================================================================
# 3. MultiSpectrum — all columns, no technique
# =============================================================================

ms = Spectrum.read(data_file, reader="spec")
print("=== MultiSpectrum (no technique) ===")
print(f"  {len(ms.fields)} fields: {[f.name for f in ms.fields]}")
print()

# =============================================================================
# 4. Figure: show all reader capabilities
# =============================================================================

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle("Fe foil EXAFS — SpecDATReader demo", fontsize=14)

# --- (0,0) Raw signals from MultiSpectrum, scan 1 ---
ax = axes[0, 0]
ax2 = ax.twinx()
ms[I0_col].plot(ax=ax, color="C0", label=f"I0 ({I0_col})")
ms[It_col].plot(ax=ax2, color="C1", label=f"It ({It_col})")
ax.axvline(7112, color="gray", lw=0.8, ls="--")
ax.set_title("Raw signals (scan 1)")
ax.set_ylabel(f"I0 counts [{I0_col}]", color="C0")
ax2.set_ylabel(f"It counts [{It_col}]", color="C1")
ax.tick_params(axis="y", colors="C0")
ax2.tick_params(axis="y", colors="C1")
ax.legend(loc="upper left", fontsize=8)
ax2.legend(loc="lower right", fontsize=8)

# --- (0,1) Scan selection: read each scan individually ---
ax = axes[0, 1]
scans = {}
for n in (1, 2, 3):
    scans[n] = Spectrum.read(
        data_file,
        reader="spec",
        technique="XAS",
        y_name=It_col,
        ref_name=I0_col,
        scan_numbers=n,
    )
    scans[n].plot(ax=ax, label=f"scan {n}", alpha=0.75)
ax.axvline(7112, color="gray", lw=0.8, ls="--")
ax.set_title("Individual scans (scan_numbers=n)")
ax.legend(fontsize=8)

# --- (0,2) Averaging: single vs averaged ---
ax = axes[0, 2]
xas_avg = Spectrum.read(
    data_file,
    reader="spec",
    technique="XAS",
    y_name=It_col,
    ref_name=I0_col,
    average_scans=True,
)
scans[1].plot(ax=ax, color="C0", alpha=0.5, label="scan 1 only")
xas_avg.plot(
    ax=ax, color="C2", label=f"average of {len(xas_avg.metadata['scan_numbers'])} scans"
)
ax.axvline(7112, color="gray", lw=0.8, ls="--")
ax.set_title("Single vs averaged (average_scans=True)")
ax.legend(fontsize=8)

# --- (1,0) Full energy range, averaged ---
ax = axes[1, 0]
xas_avg.plot(ax=ax, color="C2")
ax.axvline(7112, color="gray", lw=0.8, ls="--", label="Fe K-edge (7112 eV)")
ax.set_title("Full EXAFS range (averaged)")
ax.legend(fontsize=8)

# --- (1,1) Edge region zoomed ---
ax = axes[1, 1]
xas_avg.plot(ax=ax, color="C2")
ax.axvline(7112, color="gray", lw=0.8, ls="--", label="Fe K-edge (7112 eV)")
ax.set_xlim(7050, 7400)
e, mu = xas_avg.x, xas_avg.y
mask = (e > 7050) & (e < 7400)
ax.set_ylim(mu[mask].min() - 0.02, mu[mask].max() + 0.02)
ax.set_title("Edge region zoom (averaged)")
ax.legend(fontsize=8)

# --- (1,2) Metadata summary as text panel ---
ax = axes[1, 2]
ax.axis("off")
lines = [
    f"name:     {xas_avg.name}",
    f"tstamp:   {datetime.datetime.fromtimestamp(xas_avg.tstamp):%Y-%m-%d %H:%M:%S}",
    f"duration: {xas_avg.duration:.0f} s / scan",
    f"reader:   {type(xas_avg.reader).__name__}",
    "",
    f"scan_numbers: {xas_avg.metadata['scan_numbers']}",
    "scan_command:",
    f"  {xas_avg.metadata['scan_command']}",
    f"count_time_ms: {xas_avg.metadata['count_time_ms']}",
    "",
    "beamline_file:",
    f"  {xas_avg.metadata['beamline_file']}",
    "",
    f"comment: {xas_avg.metadata['comment']}",
]
ax.text(
    0.05,
    0.95,
    "\n".join(lines),
    transform=ax.transAxes,
    fontsize=8,
    verticalalignment="top",
    fontfamily="monospace",
    bbox=dict(boxstyle="round", facecolor="whitesmoke", alpha=0.8),
)
ax.set_title("Stored metadata")

fig.tight_layout()
plt.show()
