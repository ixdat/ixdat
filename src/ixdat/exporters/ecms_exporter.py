from .ec_exporter import ECExporter
from ..tools import deprecate


class ECMSExporter(ECExporter):
    """A CSVExporter that by default exports potential, current, selector, and all MID"""

    @property
    def default_export_columns(self):
        """The default v_list for ECExporter is V_str, J_str, and sel_str"""
        v_list = (
            ECExporter(measurement=self.measurement).default_export_columns
            + self.measurement.mass_list
        )

        return v_list

    @deprecate("0.2.0", "use `columns` instead", "0.3.0", kwarg_name="v_list")
    def export(
        self,
        path_to_file=None,
        measurement=None,
        v_list=None,
        columns=None,
        tspan=None,
        mass_list=None,
        mol_list=None,
        time_step=None,
    ):
        """Export a given measurement to a specified file.

        To improve flexibility with inheritance, this method allocates its work to:
        - CSVExporter.prepare_header_and_data()
        - CSVExporter.write_header()
        - CSVExporter.write_data()

        Args:
            measurement (Measurement): The measurement to export.
                Defaults to self.measurement.
                TODO: remove this kwarg. See conversation here:
                   https://github.com/ixdat/ixdat/pull/30/files#r810926968
            path_to_file (Path): The path to the file to write. If it has no suffix,
                a .csv suffix is appended. Defaults to f"{measurement.name}.csv"
            columns (list of str): The names of the data series to include. Defaults in
                CSVExporter to all VSeries and TSeries in the measurement. This default
                may be overwritten in inheriting exporters.
            v_list: DEPRECATED. Use `columns` instead.
            tspan (timespan): The timespan to include in the file, defaults to all of it
            time_step (float): Optional. The time spacing between data points. Can be
                used to reduce file size.
            mass_list (list of str): Names of masses to export. Defaults to all.
            mol_list (list of str or `ECMSCalResult`): Names of mols for which to export
                quantified data.
        """
        columns = columns or v_list  # deal with deprecated argument
        if not columns:
            if mass_list:
                columns = ECExporter(measurement=self.measurement).default_export_columns
            else:
                columns = self.default_export_columns
        if mass_list:
            columns += mass_list
        if mol_list:
            columns += [f"n_dot_{mol}" for mol in mol_list]
        return super().export(
            path_to_file, measurement, columns, tspan, time_step=time_step
        )
