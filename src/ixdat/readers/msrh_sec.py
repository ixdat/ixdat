from pathlib import Path  # noqa
import numpy as np
import pandas as pd
from .reading_tools import prompt_for_tstamp
from ..techniques import TECHNIQUE_CLASSES
from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..techniques.analysis_tools import calc_t_using_scan_rate


class MsrhSECReader:
    """A reader for SEC saved in three files: spectra vs v; wavelengths; current vs v"""

    def read(
        self,
        path_to_file,
        path_to_ref_spec_file,
        path_to_V_J_file,
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
            path_to_V_J_file (Path or str): The full path to the file containing the
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
        jv_df = pd.read_csv(path_to_V_J_file, names=["v", "j"])

        # The spectra need (i) the first colum with the arbitrary counter to be
        #    discarded and (ii) axes switched so that wavelength is the inner axis
        #    (axis=1). The outer axis (axis=0) then spans potential or, eq., time:
        spectra = sec_df.to_numpy()[:, 1:].swapaxes(0, 1)
        # The potential comes from the keys of that data, discarding the first column:
        v = np.array([float(key) for key in sec_df.keys()])[1:]
        # We get time from this potential and the scan rate, with a helper function:
        t = calc_t_using_scan_rate(v, dvdt=scan_rate * 1e-3)
        # If they didn't provide a tstamp, we have to prompt for it.
        tstamp = tstamp or prompt_for_tstamp(path_to_file)
        # Ready to define the measurement's TimeSeries:
        tseries = TimeSeries(
            "time from scan rate", unit_name="s", data=t, tstamp=tstamp
        )

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
        # The spectra span a space defined by time and wavelength:
        spectra = Field(
            name="spectra",
            unit_name="counts",
            axes_series=[tseries, wl_series],
            data=spectra,
        )

        # Now we process the current and potential:
        v_0 = jv_df["v"].to_numpy()  # ... but we'll actually use v from the sec data
        j_0 = jv_df["j"].to_numpy() * 1e3  # 1e3 converts [A] to [mA]
        if v_0[0] > v_0[-1]:  # Need the potential in the EC file to be increasing:
            v_0 = np.flip(v_0)
            j_0 = np.flip(j_0)
        # Since the "real" potential is in the sec data, we need to interpolate the
        #   current onto it:
        j = np.interp(v, v_0, j_0)
        # and now we're ready to define the electrochemical DataSeries:
        v_series = ValueSeries("raw potential / [V]", "V", v, tseries=tseries)
        j_series = ValueSeries("raw current / [mA]", "mA", j, tseries=tseries)

        # put all our DataSeries together:
        series_list = [
            tseries,
            v_series,
            j_series,
            wl_series,
            reference,
            spectra,
        ]

        # Figure out which measurement class to return. Use S-EC unless this read
        #   function is provided an even more specific technique class:
        measurement_class = TECHNIQUE_CLASSES["S-EC"]
        if issubclass(cls, measurement_class):
            measurement_class = cls

        # and initiate the measurement:
        measurement = measurement_class(
            name=str(path_to_file),
            tstamp=tstamp,
            series_list=series_list,
            raw_potential_names=(v_series.name,),
            raw_current_names=(j_series.name,),
        )

        return measurement


class MsrhSECDecayReader:
    def read(
        self,
        path_to_file,
        path_to_ref_spec_file,
        path_to_t_J_file,
        path_to_t_V_file,
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
            path_to_t_V_file (Path or str): The full path to the file containing the
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
        t_V_df = pd.read_csv(path_to_t_V_file, names=["t", "V"])
        t_J_df = pd.read_csv(path_to_t_J_file, names=["t", "J"])

        t_and_spectra = sec_df.to_numpy()
        spectra = t_and_spectra[:, 1:].swapaxes(0, 1)
        t_spectra = np.array([float(key) for key in sec_df.keys()])[1:]

        wl = ref_df["wavelength"].to_numpy()
        excess_wl_points = len(wl) - spectra.shape[1]
        wl = wl[excess_wl_points:]
        ref_signal = ref_df["counts"].to_numpy()[excess_wl_points:]

        wl_series = DataSeries("wavelength / [nm]", "nm", wl)
        reference = Field(
            name="reference",
            unit_name="counts",
            axes_series=[wl_series],
            data=np.array(ref_signal),
        )

        v = t_V_df["V"].to_numpy()
        t_v = t_V_df["t"].to_numpy()
        j = t_J_df["J"].to_numpy() * 1e3  # Convert [A] to [mA]
        t_j = t_J_df["t"].to_numpy()

        tstamp = tstamp or prompt_for_tstamp(path_to_file)

        tseries_j = TimeSeries("t for current", "s", data=t_j, tstamp=tstamp)
        tseries_v = TimeSeries("t for potential", "s", data=t_v, tstamp=tstamp)
        tseries_spectra = TimeSeries("t for spectra", "s", t_spectra, tstamp)
        v_series = ValueSeries("raw potential / [V]", "V", v, tseries=tseries_v)
        j_series = ValueSeries("raw current / [mA]", "mA", j, tseries=tseries_j)
        spectra = Field(
            name="spectra",
            unit_name="counts",
            axes_series=[tseries_spectra, wl_series],
            data=spectra,
        )
        series_list = [
            tseries_j,
            tseries_v,
            tseries_spectra,
            v_series,
            j_series,
            wl_series,
            reference,
            spectra,
        ]

        measurement_class = TECHNIQUE_CLASSES["S-EC"]
        if issubclass(cls, measurement_class):
            measurement_class = cls

        measurement = measurement_class(
            name=str(path_to_file),
            tstamp=tstamp,
            series_list=series_list,
            raw_potential_names=(v_series.name,),
            raw_current_names=(j_series.name,),
        )

        return measurement
