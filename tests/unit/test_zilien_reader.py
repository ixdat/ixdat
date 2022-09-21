from ixdat import Measurement


DATA_DIR = Path(__file__).parent.parent.parent / "submodules" / "ixdat-large-test-files"


class TestZilienTSVReader:
    def test_first(self):
        m = Measurement.read(
            DATA_DIR
            / "zilien_version_2/2022-09-19 11_27_59 test-1/2022-09-19 11_27_59 test-21.tsv",
            reader="zilien",
        )
        print(m.series_list)
