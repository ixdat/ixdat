from pathlib import Path
import numpy as np
import pandas as pd
from .reading_tools import prompt_for_tstamp
from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..techniques.analysis_tools import calc_t_using_scan_rate


class MsrhSECReader:
    def read(
        self,
        path_to_file,
        path_to_wl_file,
        path_to_jv_file,
        scan_rate,
        tstamp=None,
        cls=None,
    ):
        if not cls:
            from ..techniques.spectroelectrochemistry import SpectroECMeasurement

            cls = SpectroECMeasurement

        path_to_file = Path(path_to_file)
        path_to_wl_file = Path(path_to_wl_file)
        path_to_jv_file = Path(path_to_jv_file)

        sec_df = pd.read_csv(path_to_file)
        whitelight_df = pd.read_csv(path_to_wl_file, names=["wavelength", "intensity"])
        jv_df = pd.read_csv(path_to_jv_file, names=["v", "j"])

        spectra = sec_df.to_numpy()[:, 1:].swapaxes(0, 1)

        wl = whitelight_df["wavelength"].to_numpy()
        excess_wl_points = len(wl) - spectra.shape[1]
        wl = wl[excess_wl_points:]
        whitelight = whitelight_df["intensity"].to_numpy()[excess_wl_points:]

        wl_series = DataSeries("wavelength", "nm", wl)
        white_series = Field(
            "white_light",
            "counts",
            axes_series=[wl_series],
            data=np.array([whitelight]),
        )

        v = jv_df["v"].to_numpy()
        j = jv_df["j"].to_numpy()
        excess_jv_points = len(v) - spectra.shape[0]
        v = v[:-excess_jv_points]
        j = j[:-excess_jv_points]
        t = calc_t_using_scan_rate(v, dvdt=scan_rate * 1e-3)

        tstamp = tstamp or prompt_for_tstamp(path_to_file)
        tseries = TimeSeries(
            "time from scan rate", unit_name="s", data=t, tstamp=tstamp
        )
        v_series = ValueSeries("raw_potential", "V", v, tseries=tseries)
        j_series = ValueSeries("raw_current", "mA", j, tseries=tseries)
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
            white_series,
            spectra,
        ]

        measurement = cls(
            name=str(path_to_file),
            tstamp=tstamp,
            series_list=series_list,
            raw_potential_names=(v_series.name,),
            raw_current_names=(j_series.name,),
        )

        return measurement
