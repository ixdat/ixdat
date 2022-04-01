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
        aliases = super().aliases.copy()
        for name, prop_names in [
            ("t", ("t_name",)),
            ("raw_potential", ("E_name", "U_name")),
            ("raw_current", ("I_name", "J_name")),
            ("selector", ("selector_name",)),
        ]:
            # This is bit complex because the essential series for an ECMeasureemnt
            #   are t, raw_current, and raw_potential, but by default the calibrated
            #   potential and normalized current are exported if available. So
            #   we need to go through and make sure the reader of the exported .csv can
            #   find a "raw_current" and "raw_potential" using aliases, even if it
            #   isn't actually "raw". And at the same time we have to avoid circular
            #   lookups in aliases. Here goes:
            for prop_name in prop_names:
                name_in_measurement = getattr(self.measurement, prop_name)
                if (
                    name_in_measurement in self.columns
                    and (name not in aliases or name_in_measurement not in aliases[name])
                    and name_in_measurement != name
                ):
                    aliases[name] = (name_in_measurement,)
                    break
        return aliases
