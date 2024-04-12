from pathlib import Path  # noqa
import numpy as np
import pandas as pd
import time
from .reading_tools import prompt_for_tstamp
from ..techniques import TECHNIQUE_CLASSES
from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..spectra import Spectrum, SpectrumSeries
from ..techniques.analysis_tools import calc_t_using_scan_rate


class MsrhSECReader:
    """A reader for SEC saved in three files: spectra vs U; wavelengths; current vs U"""

    def read(
        self,
        path_to_file,
        path_to_ref_spec_file,
        path_to_U_J_file,
        scan_rate,
        tstamp=None,
        cls=None,
    ):
        """Read potential-dep. SEC data from 3 csv's to return a SpectroECMeasurement

        The function is well-commented so take a look at the source

        Args:
            path_to_file (Path or str): The full path to the file containing the
                spectra data. This file has voltage in the first row, and a first
                column with an arbitrary counter which has to be replaced by wavelength.
            path_to_ref_spec_file (Path or str): The full path to the file containing
                the wavelenth data, together usually with the adsorption-free spectrum.
                The length of the columns should be the same as in the spectrum data
                but in practice is a few points longer. The excess points at the starts
                of the columns are discarded.
            path_to_U_J_file (Path or str): The full path to the file containing the
                current data vs potential. The columns may be reversed in order. In the
                end the potential in the spectra file will be retained and the potential
                here used to interpolate the current onto the spectra file's potential.
            scan_rate (float): Scan rate in [mV/s]. This is used to figure out the
                measurement's time variable, as time is bizarrely not included in any
                of the data files.
            tstamp (float): Timestamp. If None, the user will be prompted for the
                measurement start time or whether to use the file creation time. This is
                necessary because tstamp is also not included in any of the files but is
                central to how ixdat organizes data. If you're sure that tstamp doesn't
                matter for you, put e.g. tstamp=1 to suppress the prompt.
            cls (Measurement subclass): The class of measurement to return. Defaults to
                SpectroECMeasurement.
        """
        # us pandas to load the data from the csv files into dataframes:
        sec_df = pd.read_csv(path_to_file)
        # ^ Note the first row, containing potential, will become the keys. The first
        #   column, containing an arbitrary counter, is included in the data.
        ref_df = pd.read_csv(path_to_ref_spec_file, names=["wavelength", "counts"])
        EI_df = pd.read_csv(path_to_U_J_file, names=["U", "J"])

        # The spectra need (i) the first colum with the arbitrary counter to be
        #    discarded and (ii) axes switched so that wavelength is the inner axis
        #    (axis=1). The outer axis (axis=0) then spans potential or, eq., time:
        spectra = sec_df.to_numpy()[:, 1:].swapaxes(0, 1)
        # The potential comes from the keys of that data, discarding the first column:
        U = np.array([float(key) for key in sec_df.keys()])[1:]
        # We get time from this potential and the scan rate, with a helper function:
        t = calc_t_using_scan_rate(U, dvdt=scan_rate * 1e-3)
        # If they didn't provide a tstamp, we have to prompt for it.
        tstamp = tstamp or prompt_for_tstamp(path_to_file)
        if tstamp == "now":
            tstamp = time.time()
        # Ready to define the measurement's TimeSeries:
        tseries = TimeSeries("time from scan rate", unit_name="s", data=t, tstamp=tstamp)

        # The wavelength comes from the reference spectrum file.
        wl = ref_df["wavelength"].to_numpy()
        excess_wl_points = len(wl) - spectra.shape[1]
        # ^ This is how many points to discard to line up with sec data
        #   (1 or 2 points in the example data).
        wl = wl[excess_wl_points:]
        ref_signal = ref_df["counts"].to_numpy()[excess_wl_points:]

        # Now we're ready to define all the spectrum DataSeries:

        # wavelength is independent variable --> simple DataSeries
        wl_series = DataSeries("wavelength / [nm]", "nm", wl)
        # The reference spectrum spans a space defined by wavelength:
        reference = Field(
            name="reference",
            unit_name="counts",
            axes_series=[wl_series],
            data=ref_signal,
        )
        reference_spectrum = Spectrum.from_field(reference)
        # The spectra span a space defined by time and wavelength:
        spectra = Field(
            name="spectra",
            unit_name="counts",
            axes_series=[tseries, wl_series],
            data=spectra,
        )
        spectrum_series = SpectrumSeries.from_field(
            spectra, tstamp=tstamp, continuous=True
        )

        # Now we process the current and potential:
        U_0 = EI_df["U"].to_numpy()  # ... but we'll actually use U from the sec data
        J_0 = EI_df["J"].to_numpy() * 1e3  # 1e3 converts [A] to [mA]
        if U_0[0] > U_0[-1]:  # Need the potential in the EC file to be increasing:
            U_0 = np.flip(U_0)
            J_0 = np.flip(J_0)
        # Since the "real" potential is in the sec data, we need to interpolate the
        #   current onto it:
        J = np.interp(U, U_0, J_0)
        # and now we're ready to define the electrochemical DataSeries:
        U_series = ValueSeries("raw potential / [V]", "V", U, tseries=tseries)
        J_series = ValueSeries("raw current / [mA]", "mA", J, tseries=tseries)

        # put all our DataSeries together:
        series_list = [
            tseries,
            U_series,
            J_series,
        ]

        # Figure out which measurement class to return. Use SEC unless this read
        #   function is provided an even more specific technique class:
        measurement_class = TECHNIQUE_CLASSES["EC-Optical"]
        if issubclass(cls, measurement_class):
            measurement_class = cls

        # and initiate the measurement:
        measurement = measurement_class(
            name=str(path_to_file),
            tstamp=tstamp,
            series_list=series_list,
            aliases={
                "raw_potential": (U_series.name,),
                "raw_current": (J_series.name,),
                "t": (tseries.name,),
            },
            spectrum_series=spectrum_series,
            reference_spectrum=reference_spectrum,
        )

        return measurement


class MsrhSECDecayReader:
    def read(
        self,
        path_to_file,
        path_to_ref_spec_file,
        path_to_t_J_file,
        path_to_t_U_file,
        tstamp=None,
        cls=None,
    ):
        """Read time-dependent SEC data from 4 csv's to return a SpectroECMeasurement

        The function works in a very similar way to MsrhSECReader.read().

        Args:
            path_to_file (Path or str): The full path to the file containing the
                spectra data. This file has time in the first row, and a first
                column with an arbitrary counter which has to be replaced by wavelength.
            path_to_ref_spec_file (Path or str): The full path to the file containing
                the wavelenth data, together usually with the adsorption-free spectrum.
                The length of the columns should be the same as in the spectrum data
                but in practice is a few points longer. The excess points at the starts
                of the columns are discarded.
            path_to_t_U_file (Path or str): The full path to the file containing the
                potential data vs time.
            path_to_t_J_file (Path or str): The full path to the file containing the
                current data vs time.
            tstamp (float): Timestamp. If None, the user will be prompted for the
                measurement start time or whether to use the file creation time. This is
                necessary because tstamp is also not included in any of the files but is
                central to how ixdat organizes data. If you're sure that tstamp doesn't
                matter for you, put e.g. tstamp=1 to suppress the prompt.
            cls (Measurement subclass): The class of measurement to return. Defaults to
                SpectroECMeasurement.
        """

        sec_df = pd.read_csv(path_to_file)
        ref_df = pd.read_csv(path_to_ref_spec_file, names=["wavelength", "counts"])
        t_U_df = pd.read_csv(path_to_t_U_file, names=["t", "U"])
        t_J_df = pd.read_csv(path_to_t_J_file, names=["t", "J"])

        t_and_spectra = sec_df.to_numpy()
        spectra = t_and_spectra[:, 1:].swapaxes(0, 1)
        t_spectra = np.array([float(key) for key in sec_df.keys()])[1:]

        wl = ref_df["wavelength"].to_numpy()
        excess_wl_points = len(wl) - spectra.shape[1]
        wl = wl[excess_wl_points:]
        ref_signal = ref_df["counts"].to_numpy()[excess_wl_points:]

        wl_series = DataSeries("wavelength / [nm]", "nm", wl)
        reference_spectrum = Spectrum.from_field(
            Field(
                name="reference",
                unit_name="counts",
                axes_series=[wl_series],
                data=np.array(ref_signal),
            ),
            tstamp=tstamp,
        )

        U = t_U_df["U"].to_numpy()
        t_U = t_U_df["t"].to_numpy()
        J = t_J_df["J"].to_numpy() * 1e3  # Convert [A] to [mA]
        t_E = t_J_df["t"].to_numpy()

        tstamp = tstamp or prompt_for_tstamp(path_to_file)
        if tstamp == "now":
            tstamp = time.time()

        tseries_J = TimeSeries("t for current", "s", data=t_E, tstamp=tstamp)
        tseries_U = TimeSeries("t for potential", "s", data=t_U, tstamp=tstamp)
        tseries_spectra = TimeSeries("t for spectra", "s", t_spectra, tstamp)
        U_series = ValueSeries("raw potential / [V]", "V", U, tseries=tseries_U)
        J_series = ValueSeries("raw current / [mA]", "mA", J, tseries=tseries_J)
        spectrum_series = SpectrumSeries.from_field(
            Field(
                name="spectra",
                unit_name="counts",
                axes_series=[tseries_spectra, wl_series],
                data=spectra,
            ),
            tstamp=tstamp,
        )
        series_list = [
            tseries_J,
            tseries_U,
            tseries_spectra,
            U_series,
            J_series,
            wl_series,
        ]

        measurement_class = TECHNIQUE_CLASSES["EC-Optical"]
        if issubclass(cls, measurement_class):
            measurement_class = cls

        measurement = measurement_class(
            name=str(path_to_file),
            tstamp=tstamp,
            series_list=series_list,
            aliases={
                "raw_potential": (U_series.name,),
                "raw_current": (J_series.name,),
                "t": (tseries_U.name,),
            },
            spectrum_series=spectrum_series,
            reference_spectrum=reference_spectrum,
        )

        return measurement
