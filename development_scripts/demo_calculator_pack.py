"""
1) Read an EC measurement (Biologic .mpt)
2) Add and apply ECCalibration (compute 'potential' and 'current')
3) Pack calculators to file (exclude Indexer)
4) Read pack back
5) Read another measurement (can be same file for demo)
6) Attach pack to target and APPLY again
"""

from __future__ import annotations
from pathlib import Path

from ixdat import Measurement
from ixdat.calculators.ec_calculators import ECCalibration
from ixdat.calculators.packs import CalculatorPack

# --- CONFIG -------------------------------------------------------------------
TEST_FILE = Path(__file__).parent / "../test_data/biologic/Pt_poly_cv.mpt"
PACK_PATH = Path(__file__).parent / "demo_pack.ixpack.json"


def compute_ec_series(meas: Measurement):
    """Trigger calculation of EC-calibrated series to prove calibration is applied."""
    pot = meas["potential"].data
    cur = meas["current"].data
    print(f"\tpotential: len={len(pot)}, min={pot.min():.3g}, max={pot.max():.3g}")
    print(f"\tcurrent:   len={len(cur)}, min={cur.min():.3g}, max={cur.max():.3g}")


def main():
    # 1) Read a real EC measurement
    print("\n[1] Reading source measurement...")
    src = Measurement.read(TEST_FILE, reader="biologic")
    print(
        "\ttechnique:",
        getattr(src, "technique", None) or getattr(src, "techniques", None),
    )
    print("\tseries:", [s.name for s in src.series_list][:6], "...")
    print("\n[2] EC series before calibration...")
    compute_ec_series(src)

    # 2) Add only ECCalibration (EC-only file)
    print("\n[3] Adding ECCalibration...")
    ec_cal = ECCalibration(name="ec", A_el=0.196, RE_vs_RHE=0.2, R_Ohm=0.0)
    src.add_calculator(ec_cal)

    print("\tCalculators now on src:")
    cal_list = (
        src.calculators.values()
        if isinstance(src.calculators, dict)
        else src.calculators
    )
    for c in cal_list:
        print("\t-", c)

    # 3) APPLY EC calibration by computing series
    print("\n[4] Computing calibrated EC series on src:")
    compute_ec_series(src)

    # 4) Build a CalculatorPack and write to file (exclude Indexer by default)
    print("\n[5] Building CalculatorPack (excluding Indexer) and writing to file...")
    pack = CalculatorPack.from_measurement(
        src,
        name="demo_pack",
        notes="EC-only demo",
        exclude=["Indexer", "ixdat.calculators.indexer.Indexer"],
    )
    print(pack.summary())

    pack.to_file(str(PACK_PATH))
    print("\tWrote:", PACK_PATH.resolve())

    # 5) Read the pack back
    print("\n[6] Reading pack back from file...")
    pack2 = CalculatorPack.read(path=str(PACK_PATH))
    print(pack2.summary())

    # 6) Read target measurement (reuse same file for demo) and attach pack
    print("\n[7] Reading target measurement...")
    tgt = Measurement.read(TEST_FILE, reader="biologic")

    print("\tAttaching pack to target...")
    attached = pack2.attach_to(
        tgt,
        on_conflict="replace",
    )
    print("\tAttached calculators:")
    for c in attached:
        print("\t-", c)

    print("\n[8] Computing calibrated EC series on target:")
    compute_ec_series(tgt)

    print("\nDone")


if __name__ == "__main__":
    main()
