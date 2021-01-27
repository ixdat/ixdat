from .csv_exporter import CSVExporter


class ECExporter(CSVExporter):
    @property
    def default_v_list(self):
        return [
            # self.measurement.t_str,
            self.measurement.V_str,
            self.measurement.J_str,
            self.measurement.sel_str,
        ]
