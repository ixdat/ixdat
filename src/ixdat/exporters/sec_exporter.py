from .csv_exporter import CSVExporter
from .ec_exporter import ECExporter
from .spectrum_exporter import SpectrumExporter, SpectrumSeriesExporter


class SECExporter(CSVExporter):
    """Adds to CSVExporter the export of the Field with the SEC spectra"""

    def __init__(self, measurement, delim=",\t"):
        super().__init__(measurement, delim=delim)
        self.reference_exporter = SpectrumExporter(measurement.reference_spectrum)
        self.spectra_exporter = SpectrumSeriesExporter(measurement.spectrum_series)

    @property
    def default_v_list(self):
        """The default v_list for SECExporter is that from EC and tracked wavelengths"""
        v_list = (
            ECExporter(measurement=self.measurement).default_v_list
            + self.measurement.tracked_wavelengths
        )
        return v_list

    def prepare_header_and_data(self, measurement, v_list, tspan):
        """Add lines pointing to the 'spectra' and 'reference' spectrum"""
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
