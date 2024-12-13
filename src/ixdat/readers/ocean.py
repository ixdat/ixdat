import pandas as pd
from ..data_series import Field, DataSeries


class OceanTxTReader():
    """
    Class for reading .txt files obtained from Ocean Optics spectrometers.
    """
    def read(self, path_to_file, cls=None, **kwargs):
        df = pd.read_csv(path_to_file, sep="\t",
                         skiprows=14, names=["Wavelength", "Absorbance"])

        df["Wavelength"] = df["Wavelength"].str.replace(',', '.').astype(float)
        df["Absorbance"] = df["Absorbance"].str.replace(',', '.').astype(float)

        wl = DataSeries(name="Wavelength", unit_name="nm",
                        data=df["Wavelength"].to_numpy())
        field = Field(axes_series=[wl], name="Absorbance",
                      unit_name="a.u.", data=df["Absorbance"].to_numpy())

        return cls.from_field(field, **kwargs)
