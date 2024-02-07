import pandas as pd
from collections import OrderedDict


class SpectrumExporter:
    """An ixdat CSV exporter for spectra. Uses pandas."""

    def __init__(self, spectrum, delim=","):
        """Initiate the SpectrumExporter.

        Args:
            spectrum (Spectrum): The spectrum to export by default
            delim (char): The separator for the .csv file. Note that this cannot be
                the ",\t" used by ixdat's main exporter since pandas only accepts single
                character delimiters.
        """
        self.spectrum = spectrum
        self.delim = delim

    def export(self, path_to_file, spectrum=None):
        """Export spectrum to path_to_file.

        Args:
            path_to_file (str or Path): The path of the file to export to. Note that if a
                file already exists with this path, it will be overwritten.
            spectrum (Spectrum): The spectrum to export if different from self.spectrum
                TODO: remove this kwarg. See conversation here:
                   https://github.com/ixdat/ixdat/pull/30/files#r810926968
        """
        spectrum = spectrum or self.spectrum
        df = pd.DataFrame({spectrum.x_name: spectrum.x, spectrum.y_name: spectrum.y})

        header_lines = []
        for attr in ["name", "technique", "tstamp", "backend_name", "id"]:
            line = f"{attr} = {getattr(spectrum, attr)}\n"
            header_lines.append(line)

        header_lines.append("\n")

        # Insert a line, after the first two lines, saying how long the header is.
        N_header_lines = len(header_lines) + 2
        # The `+ 2` is for the header length line and the column header line.
        header_length_line = f"N_header_lines = {N_header_lines}\n"
        header_lines = header_lines[:2] + [header_length_line] + header_lines[2:]

        with open(path_to_file, "w") as f:
            f.writelines(header_lines)
        with open(path_to_file, "a") as f:
            df.to_csv(f, index=False, sep=self.delim, lineterminator="\n")

        print(f"wrote {path_to_file}!")


class SpectrumSeriesExporter:
    """An exporter for ixdat spectrum series."""

    def __init__(self, spectrum_series, delim=","):
        """Initiate the SpectrumSeriesExporter.

        Args:
            spectrum_series (SpectrumSeries): The spectrum to export by default
            delim (char): The separator for the .csv file. Note that this cannot be
                the ",\t" used by ixdat's main exporter since pandas only accepts single
                character delimiters.
        """
        self.spectrum_series = spectrum_series
        self.delim = delim

    def export(self, path_to_file=None, spectrum_series=None, spectra_as_rows=True):
        """Export spectrum series to path_to_file.

        Args:
            spectrum_series (Spectrum): The spectrum_series to export if different from
                self.spectrum_series
                TODO: remove this kwarg. See conversation here:
                   https://github.com/ixdat/ixdat/pull/30/files#r810926968
            path_to_file (str or Path): The path of the file to export to. Note that if a
                file already exists with this path, it will be overwritten.
            spectra_as_rows (bool): This specifies the orientation of the data exported.
                If True, the scanning variabe (e.g. wavelength) increases to the right
                and the time variable increases downward. If False, the scanning
                variable increases downwards and the time variable increases to the
                right. Either way it is clarified in the file header. Defaults to True.
        """

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

        header_lines.append("\n")

        # Insert a line, after the first two lines, saying how long the header is.
        N_header_lines = len(header_lines) + 2
        # The `+ 2` is for the header length line and the column header line.
        header_length_line = f"N_header_lines = {N_header_lines}\n"
        header_lines = header_lines[:2] + [header_length_line] + header_lines[2:]

        with open(path_to_file, "w") as f:
            f.writelines(header_lines)
        with open(path_to_file, "a") as f:
            df.to_csv(f, index=False, sep=self.delim, lineterminator="\n")

        print(f"wrote {path_to_file}!")
