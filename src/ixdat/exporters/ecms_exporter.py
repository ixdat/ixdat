from .csv_exporter import CSVExporter
from .ec_exporter import ECExporter


class ECMSExporter(CSVExporter):
    """A CSVExporter that by default exports potential, current, selector, and all MID"""

    @property
    def default_v_list(self):
        """The default v_list for ECExporter is V_str, J_str, and sel_str"""
        v_list = (
            ECExporter(measurement=self.measurement).default_v_list
            + self.measurement.mass_list
        )

        return v_list
