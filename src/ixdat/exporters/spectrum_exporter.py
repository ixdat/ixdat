import pandas as pd
from collections import OrderedDict


class SpectrumExporter:
    def __init__(self, spectrum, delim=","):
        self.spectrum = spectrum
        self.delim = delim

    def export(self, spectrum, path_to_file):
        spectrum = spectrum or self.spectrum
        df = pd.DataFrame({spectrum.x_name: spectrum.x, spectrum.y_name: spectrum.y})

        header_lines = []
        for attr in ["name", "technique", "tstamp", "backend_name", "id"]:
            line = f"{attr} = {getattr(spectrum, attr)}\n"
            header_lines.append(line)

        N_header_lines = len(header_lines) + 3
        header_lines.append(f"N_header_lines = {N_header_lines}\n")
        header_lines.append("\n")
        # header_lines.append("".join([(key + self.delim) for key in df.keys()]))
        df.to_csv(path_to_file, index=False, sep=self.delim)

        with open(path_to_file, "w") as f:
            f.writelines(header_lines)
        with open(path_to_file, "a") as f:
            df.to_csv(f, index=False, sep=self.delim, line_terminator="\n")

        print(f"wrote {path_to_file}!")


class SpectrumSeriesExporter:
    def __init__(self, spectrum_series, delim=","):
        self.spectrum_series = spectrum_series
        self.delim = delim

    def export(self, spectrum_series=None, path_to_file=None, spectra_as_rows=True):

        spectrum_series = spectrum_series or self.spectrum_series

        field = spectrum_series.field
        data = field.data
        tseries, xseries = spectrum_series.field.axes_series
        t = tseries.t + tseries.tstamp - spectrum_series.tstamp
        x = xseries.data

        header_lines = []
        for attr in ["name", "technique", "tstamp", "backend_name", "id"]:
            line = f"{attr} = {getattr(spectrum_series, attr)}\n"
            header_lines.append(line)

        header_lines.append(
            f"values are y='{field.name}' with units [{field.unit_name}]\n"
        )

        if spectra_as_rows:  # columns are ValueSeries
            data_as_list_of_tuples = [(spectrum_series.t_name, t)] + [
                (x_i, data[:, i]) for i, x_i in enumerate(x)
            ]
            df = pd.DataFrame(OrderedDict(data_as_list_of_tuples))
            header_lines.append(
                f"first row is x='{xseries.name}' with units [{xseries.unit_name}]\n"
            )
            header_lines.append(
                f"first column is t='{tseries.name}' with units [{tseries.unit_name}]\n"
            )
        else:  # spectra as columns. rows are ValueSeries
            data_as_list_of_tuples = [(spectrum_series.x_name, x)] + [
                (t_i, data[i, :]) for i, t_i in enumerate(t)
            ]
            df = pd.DataFrame(OrderedDict(data_as_list_of_tuples))
            header_lines.append(
                f"first row is t='{tseries.name}' with units [{tseries.unit_name}]\n"
            )
            header_lines.append(
                f"first column is x='{xseries.name}' with units [{xseries.unit_name}]\n"
            )

        N_header_lines = len(header_lines) + 3
        header_lines.append(f"N_header_lines = {N_header_lines}\n")
        header_lines.append("\n")

        with open(path_to_file, "w") as f:
            f.writelines(header_lines)
        with open(path_to_file, "a") as f:
            df.to_csv(f, index=False, sep=self.delim, line_terminator="\n")

        print(f"wrote {path_to_file}!")
