from .ec_exporter import ECExporter
from .spectrum_exporter import SpectrumExporter, SpectrumSeriesExporter


class SECExporter(ECExporter):
    """Adds to CSVExporter the export of the Field with the SEC spectra"""

    def __init__(self, measurement, delim=","):
        super().__init__(measurement, delim=delim)
        # FIXME: The lines below don't work because this __init__ gets called before
        #   the measurement's __init__ is finished.
        # self.reference_exporter = SpectrumExporter(measurement.reference_spectrum)
        # self.spectra_exporter = SpectrumSeriesExporter(measurement.spectrum_series)
        self._reference_exporter = None
        self._spectra_exporter = None

    @property
    def reference_exporter(self):
        if not self._reference_exporter:
            self._reference_exporter = SpectrumExporter(
                self.measurement.reference_spectrum
            )
        return self._reference_exporter

    @property
    def spectra_exporter(self):
        if not self._spectra_exporter:
            self._spectra_exporter = SpectrumSeriesExporter(
                self.measurement.spectrum_series
            )
        return self._spectra_exporter

    @property
    def default_export_columns(self):
        """The default v_list for SECExporter is that from EC and tracked wavelengths"""
        columns = (
            ECExporter(measurement=self.measurement).default_export_columns
            + self.measurement.tracked_wavelengths
        )
        return columns

    def prepare_header_and_data(self, measurement, columns, tspan=None, time_step=None):
        """Do the standard ixdat csv export header preparation, plus SEC stuff.

        The SEC stuff is:
            - export the spectroelectrochemistry spectra
            - export the actual reference spectrum
            - add lines to the main file header pointing to the files with the
                above two exports.

        Args and Kwargs: see :meth:`ECExporter.prepare_header_and_data`
        """
        super().prepare_header_and_data(measurement, columns, tspan, time_step=time_step)
        path_to_spectra_file = self.path_to_file.parent / (
            self.path_to_file.stem + "_spectra.csv"
        )
        measurement = measurement or self.measurement
        self.header_lines.append(
            f"'spectrum_series' in file: '{path_to_spectra_file.name}'\n"
        )
        self.spectra_exporter.export(path_to_spectra_file)
        path_to_reference_spectrum_file = self.path_to_file.parent / (
            self.path_to_file.stem + "_reference.csv"
        )
        self.header_lines.append(
            f"'reference_spectrum' in file: '{path_to_reference_spectrum_file.name}'\n"
        )
        self.reference_exporter.export(path_to_reference_spectrum_file)

        print(f"writing {self.path_to_file}!")
