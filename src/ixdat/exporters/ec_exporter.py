from .csv_exporter import CSVExporter


class ECExporter(CSVExporter):
    """A CSVExporter that by default exports potential, current, and selector"""

    @property
    def default_export_columns(self):
        """The default v_list for ECExporter is V_str, J_str, and sel_str"""
        return [
            # self.measurement.t_name,
            self.measurement.v_name,
            self.measurement.j_name,
            self.measurement.selector_name,
        ]

    @property
    def aliases(self):
        return {
            "t": (self.measurement.t_name,),
            "raw_potential": (self.measurement.v_name,),
            "raw_current": (self.measurement.j_name,),
            "selector": (self.measurement.selector_name,),
        }
