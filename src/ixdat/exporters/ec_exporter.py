from .csv_exporter import CSVExporter


class ECExporter(CSVExporter):
    """A CSVExporter that by default exports potential, current, and selector"""

    @property
    def default_v_list(self):
        """The default v_list for ECExporter is V_str, J_str, and sel_str"""
        return [
            # self.measurement.t_str,
            self.measurement.V_str,
            self.measurement.J_str,
            self.measurement.selector_name,
        ]
