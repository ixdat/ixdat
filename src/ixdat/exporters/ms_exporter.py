from . import CSVExporter, SpectrumSeriesExporter


class MSExporter(CSVExporter):
    """By default, everything from an MS Measurement is exported."""

    pass


class MSSpectroExporter(MSExporter):
    def __init__(self, measurement, delim=","):
        super().__init__(measurement, delim=delim)
        # self.spectra_exporter = SpectrumSeriesExporter(measurement.spectrum_series)
        # FIXME: Have to do a property because this __int__ gets called before the
        #    measurement's __init__ is finished...
        self._spectra_exporter = None

    @property
    def spectra_exporter(self):
        if not self._spectra_exporter:
            self._spectra_exporter = SpectrumSeriesExporter(
                self.measurement.spectrum_series
            )
        return self._spectra_exporter

    def prepare_header_and_data(self, measurement, columns, tspan=None, time_step=None):
        """Do the standard ixdat csv export header preparation, plus SEC stuff.

        The MS Spectra stuff is:
            - export the MSSpectrumSeries
            - add a line to the main file header pointing to the spectra file

        Args and Kwargs: see :meth:`ECExporter.prepare_header_and_data`
        """
        super().prepare_header_and_data(measurement, columns, tspan, time_step=time_step)
        path_to_spectra_file = self.path_to_file.parent / (
            self.path_to_file.stem + "_spectra.csv"
        )
        self.header_lines.append(
            f"'spectrum_series' in file: '{path_to_spectra_file.name}'\n"
        )
        self.spectra_exporter.export(path_to_spectra_file)

        print(f"writing {self.path_to_file}!")
