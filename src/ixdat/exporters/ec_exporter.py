from .csv_exporter import CSVExporter


class ECExporter(CSVExporter):
    """A CSVExporter that by default exports potential, current, and selector"""

    @property
    def default_export_columns(self):
        """The default v_list for ECExporter is V_str, J_str, and sel_str"""
        return [
            # self.measurement.t_name,  # gets included automatically.
            self.measurement.U_name,
            self.measurement.J_name,
            self.measurement.selector_name,
        ]

    @property
    def aliases(self):  # TODO: Figure out if this can be deleted.
        """Ensure that the essential series are accessible as aliases."""
        aliases = self.measurement.aliases.copy()
        for name, prop_name in [
            ("t", "t_name"),
            ("raw_potential", "E_name"),
            ("raw_current", "I_name"),
            ("selector", "selector_name"),
        ]:
            name_in_measurement = getattr(self.measurement, prop_name)
            if name not in aliases and name_in_measurement != name:
                aliases[name] = (name_in_measurement,)
        return aliases
