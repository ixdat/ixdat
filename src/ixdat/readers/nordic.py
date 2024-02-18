from ..exceptions import ReadError
from ..techniques import ECMeasurement

class NordicTDMSReader:

    def __init__(self):
        self.tdms_file = None
    def read(self, path_to_file, *, cls=ECMeasurement, **kwargs):

        try:
            from nptdms import TdmsFile
        except ImportError as e:
            raise ReadError(
                "To read Nordic Electrochemistry .tdms files, ixdat uses "
                "the npTDMS package. Please install npTDMS and try again.\n"
                "see: https://nptdms.readthedocs.io/en/stable/quickstart.html "
                f"\noriginal error: \n{e}"
            )

        self.tdms_file = TdmsFile.read(path_to_file)
        name = self.tdms_file.properties["name"]
        tstamp = self.tdms_file.properties["dateTime"].astype('datetime64[s]').astype('int')


