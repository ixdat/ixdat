from .csv_exporter import CSVExporter
from .ec_exporter import ECExporter
from .spectrum_exporter import SpectrumExporter, SpectrumSeriesExporter


class SECExporter(CSVExporter):
    """Adds to CSVExporter the export of the Field with the SEC spectra"""

    def __init__(self, measurement, delim=",\t"):
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
        v_list = (
            ECExporter(measurement=self.measurement).default_v_list
            + self.measurement.tracked_wavelengths
        )
        return v_list

    aliases = ECExporter.aliases

    def prepare_header_and_data(self, measurement, v_list, tspan):
        """Do the standard ixdat csv export header preparation, plus SEC stuff.

        The SEC stuff is:
            - export the spectroelectrochemistry spectra
            - export the actual reference spectrum
            - add lines to the main file header pointing to the files with the
                above two exports.
        """
        super().prepare_header_and_data(measurement, v_list, tspan)
        path_to_spectra_file = self.path_to_file.parent / (
            self.path_to_file.stem + "_spectra.csv"
        )
        measurement = measurement or self.measurement
        self.header_lines.append(f"'spectra' in file: '{path_to_spectra_file.name}'\n")
        self.spectra_exporter.export(measurement.spectrum_series, path_to_spectra_file)
        path_to_reference_spectrum_file = self.path_to_file.parent / (
            self.path_to_file.stem + "_reference.csv"
        )
        self.header_lines.append(
            f"'reference' in file: '{path_to_reference_spectrum_file.name}'\n"
        )
        self.reference_exporter.export(
            measurement.reference_spectrum, path_to_reference_spectrum_file
        )

        print(f"writing {self.path_to_file}!")
