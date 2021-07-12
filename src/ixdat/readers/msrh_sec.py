from pathlib import Path
import numpy as np
import pandas as pd
from .reading_tools import prompt_for_tstamp
from ..techniques import TECHNIQUE_CLASSES
from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..techniques.analysis_tools import calc_t_using_scan_rate


class MsrhSECReader:
    def read(
        self,
        path_to_file,
        path_to_ref_spec_file,
        path_to_V_J_file,
        scan_rate,
        tstamp=None,
        cls=None,
    ):

        measurement_class = TECHNIQUE_CLASSES["S-EC"]
        if issubclass(cls, measurement_class):
            measurement_class = cls

        path_to_file = Path(path_to_file)
        path_to_ref_spec_file = Path(path_to_ref_spec_file)
        path_to_V_J_file = Path(path_to_V_J_file)

        sec_df = pd.read_csv(path_to_file)
        ref_df = pd.read_csv(path_to_ref_spec_file, names=["wavelength", "counts"])
        jv_df = pd.read_csv(path_to_V_J_file, names=["v", "j"])

        spectra = sec_df.to_numpy()[:, 1:].swapaxes(0, 1)

        wl = ref_df["wavelength"].to_numpy()
        excess_wl_points = len(wl) - spectra.shape[1]
        wl = wl[excess_wl_points:]
        ref_signal = ref_df["counts"].to_numpy()[excess_wl_points:]

        wl_series = DataSeries("wavelength / [nm]", "nm", wl)
        reference = Field(
            "reference",
            "counts",
            axes_series=[wl_series],
            data=np.array(ref_signal),
        )

        v_0 = jv_df["v"].to_numpy()
        j_0 = jv_df["j"].to_numpy()
        # v = v[:-excess_jv_points]   # WRONG!!!
        v = np.array([float(key) for key in sec_df.keys()])[1:]
        j = np.flip(j_0)[-len(v) :]
        t = calc_t_using_scan_rate(v, dvdt=scan_rate * 1e-3)

        tstamp = tstamp or prompt_for_tstamp(path_to_file)
        tseries = TimeSeries(
            "time from scan rate", unit_name="s", data=t, tstamp=tstamp
        )
        v_series = ValueSeries("raw potential / [V]", "V", v, tseries=tseries)
        j_series = ValueSeries("raw current / [mA]", "mA", j, tseries=tseries)
        spectra = Field(
            name="spectra",
            unit_name="counts",
            axes_series=[tseries, wl_series],
            data=spectra,
        )
        series_list = [
            tseries,
            v_series,
            j_series,
            wl_series,
            reference,
            spectra,
        ]

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
        offset=None,
        tstamp=None,
        cls=None,
    ):

        measurement_class = TECHNIQUE_CLASSES["S-EC"]
        if issubclass(cls, measurement_class):
            measurement_class = cls

        path_to_file = Path(path_to_file)
        path_to_ref_spec_file = Path(path_to_ref_spec_file)
        path_to_t_V_file = Path(path_to_t_V_file)
        path_to_t_J_file = Path(path_to_t_J_file)

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

        if False:
            # trim stuff down?
            excess_t_v_points = len(v) - spectra.shape[0]
            t_v = t_v[:-excess_t_v_points]
            v = v[:-excess_t_v_points]
            excess_t_j_points = len(j) - spectra.shape[0]
            t_j = t_j[:-excess_t_j_points]
            j = j[:-excess_t_j_points]
        if offset:
            t_j = t_j - offset  # I have no idea why but there seems a 22-second offset.

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
            v_series,
            j_series,
            wl_series,
            reference,
            spectra,
        ]

        measurement = measurement_class(
            name=str(path_to_file),
            tstamp=tstamp,
            series_list=series_list,
            raw_potential_names=(v_series.name,),
            raw_current_names=(j_series.name,),
        )

        return measurement
