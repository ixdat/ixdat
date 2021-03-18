from .csv_exporter import CSVExporter


class ECExporter(CSVExporter):
    """A CSVExporter that by default exports potential, current, and selector"""

    @property
    def default_v_list(self):
        """The default v_list for ECExporter is V_str, J_str, and sel_str"""
        v_list = [
            self.measurement.E_str,
            self.measurement.I_str,
            self.measurement.sel_str,
        ]
        if (self.measurement.RE_vs_RHE is not None) or self.measurement.R_Ohm:
            v_list.append(self.measurement.V_str)
        if self.measurement.A_el:
            v_list.append(self.measurement.J_str)
        return v_list
