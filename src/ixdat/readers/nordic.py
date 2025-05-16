from ..exceptions import ReadError
from ..techniques import ECMeasurement
from ..data_series import DataSeries, TimeSeries, ValueSeries
import numpy as np


class NordicTDMSReader:
    """Reader for the .tdms files produced by Nordic Electrochemistry potentiostat"""

    def __init__(self):
        self.tdms_file = None  # The raw TdmsFile as imported by npTDMS package

    def read(self, path_to_file, *, cls=ECMeasurement, **kwargs):
        """Read the .tdms file and return the resulting measurement object

        Requires that the npTDMS package has been installed.
        See: https://nptdms.readthedocs.io/en/stable/quickstart.html

        Args:
            path_to_file (Path or str): Path to the .tdms file to read
            cls (Measurement subclass): the class to return an object of. Defaults to
                ECMeasurement
            kwargs: Additional key-word arguments are included in the dictionary used
                to initiate the measurement by cls.from_dict().
        """
        try:
            from nptdms import TdmsFile
        except ImportError as e:
            raise ReadError(
                "To read Nordic Electrochemistry .tdms files, ixdat uses "
                "the npTDMS package. Please install npTDMS and try again.\n"
                "see: https://nptdms.readthedocs.io/en/stable/quickstart.html "
                f"\noriginal error: \n{e}"
            )

        tdms_file = TdmsFile.read(path_to_file)
        self.tdms_file = tdms_file

        name = tdms_file.properties["name"]
        tstamp = tdms_file.properties["dateTime"].astype("datetime64[s]").astype("int")

        t = tdms_file["EC"]["Time"][:]  # time in [s]
        V = tdms_file["EC"]["E"][:]  # raw potential in [V]
        i = tdms_file["EC"]["i"][:] * 1e3  # raw current in [mA]

        Z = tdms_file["EC"]["Z_E"][:]
        phase = tdms_file["EC"]["Phase_E"][:]

        tseries = TimeSeries(name="Time", unit_name="s", data=t, tstamp=tstamp)
        Vseries = ValueSeries(name="E", unit_name="V", data=V, tseries=tseries)
        Iseries = ValueSeries(name="i", unit_name="mA", data=i, tseries=tseries)

        # whoa, Z and p (impedance and phase) have different lenghts than t
        # To deal with this, we'll have to assume that they are evenly spaced in time
        # over the same timespan as potential and current:

        t_EIS = np.linspace(t[0], t[-1], num=len(Z))
        tseries_EIS = TimeSeries(
            name="Time EIS", unit_name="s", data=t_EIS, tstamp=tstamp
        )
        Zseries = ValueSeries(
            name="Z_E",
            unit_name="Ohm",
            data=Z,
            tseries=tseries_EIS,
        )
        pseries = ValueSeries(
            name="phase_E",
            unit_name=None,
            data=phase,
            tseries=tseries_EIS,
        )

        aliases = {"t": ["Time"], "raw_potential": ["E"], "raw_current": ["i"]}

        obj_as_dict = dict(
            name=name,
            technique="EC",
            tstamp=tstamp,
            series_list=[tseries, Vseries, Iseries, Zseries, pseries],
            aliases=aliases,
            reader=self,
        )
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)
