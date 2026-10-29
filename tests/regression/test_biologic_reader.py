"""Regression tests for the biologic reader"""

from pathlib import Path

from ixdat import Measurement

PATH_TO_MPT = (
    Path(__file__).parent.parent.parent / "test_data/biologic/Pt_poly_cv_CUT.mpt"
)


def test_read_set_and_cut_with_empty_mpt(tmp_path):
    """A set with a header-only .mpt file can be read and cut.

    See https://github.com/ixdat/ixdat/issues/93
    """
    lines = PATH_TO_MPT.read_text(encoding="latin-1").splitlines(keepends=True)
    n_header_lines = int(lines[1].split(":")[-1])
    (tmp_path / "set_01_CVA_C01.mpt").write_text("".join(lines), encoding="latin-1")
    (tmp_path / "set_02_CVA_C01.mpt").write_text(
        "".join(lines[:n_header_lines]), encoding="latin-1"
    )

    meas = Measurement.read_set(tmp_path / "set", suffix=".mpt", reader="biologic")
    t_start = meas.t[0]
    part = meas.cut(tspan=[t_start + 20, t_start + 30])

    assert len(part.t) > 0
