"""Classes for exporting measurement data"""
from pathlib import Path
import json
from .. import __version__


class CSVExporter:
    """The default exporter, which writes delimited measurement data row-wise to file"""

    default_export_columns = None  # Typically overwritten by inheriting Exporters
    """The names of the value series to export by default."""
    aliases = None  # This will typically be overwritten by inheriting Exporters
    """The aliases, needed for techniques with essential series that get renamed."""

    def __init__(self, measurement=None, delim=",\t"):
        """Initiate the exported with a measurement (Measurement) and delimiter (str)"""
        self.measurement = measurement
        self.delim = delim
        self.header_lines = None
        self.s_list = None
        self.columns_data = None
        self.path_to_file = None

    def export(self, path_to_file=None, measurement=None, columns=None, tspan=None):
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
            tspan (timespan): The timespan to include in the file, defaults to all of it
        """
        measurement = measurement or self.measurement
        if not path_to_file:
            path_to_file = f"{measurement.name}.csv"
        if isinstance(path_to_file, str):
            path_to_file = Path(path_to_file)
        if not path_to_file.suffix:
            path_to_file = path_to_file.with_suffix(".csv")
        self.path_to_file = path_to_file
        self.prepare_header_and_data(measurement, columns, tspan)
        self.prepare_column_header()
        self.write_header()
        self.write_data()

    def prepare_header_and_data(self, measurement, v_list, tspan):
        """Prepare self.header_lines to include metadata and value-time pairs

        Args:
            measurement (Measurement): The measurement being exported
            v_list (list of str): The names of the ValueSeries to include
            tspan (timespan): The timespan of the data to include in the export
        """
        columns_data = {}
        # list of the value names to export:
        v_list = v_list or self.default_export_columns or list(measurement.value_names)
        s_list = []  # list of the series names to export.
        # s_list will also include names of TimeSeries.

        timecols = {}  # Will be {time_name: value_names}, for the header.
        for v_name in v_list:
            # Collect data and names for each ValueSeries and TimeSeries
            t_name = measurement[v_name].tseries.name
            t, v = measurement.grab(v_name, tspan=tspan)
            if t_name in timecols:
                # We've already collected the data for this time column
                timecols[t_name].append(v_name)
            else:
                # New time column. Collect its data and add it to the timecols.
                columns_data[t_name] = t
                s_list.append(t_name)
                timecols[t_name] = [v_name]
            columns_data[v_name] = v
            s_list.append(v_name)

        header_lines = []
        ixdat_version_line = f"ixdat version = {__version__}\n"
        header_lines.append(ixdat_version_line)
        for attr in ["name", "technique", "tstamp", "backend_name", "id"]:
            line = f"{attr} = {getattr(measurement, attr)}\n"
            header_lines.append(line)
            # TODO: This should be more automated... the exporter should put all
            #    the appropriate metadata attributes of the object, read from its
            #    table definition, in the header.
        for t_name, v_names in timecols.items():
            # Header includes a line for each time column stating which values use it:
            line = (
                f"timecol '{t_name}' for: "
                + " and ".join([f"'{v_name}'" for v_name in v_names])
                + "\n"
            )
            header_lines.append(line)
        if self.aliases:
            # For now, aliases is nice after the timecol lines. But see the to-do above.
            aliases_line = f"aliases = {json.dumps(self.aliases)}\n"
            header_lines.append(aliases_line)
        self.header_lines = header_lines
        self.s_list = s_list
        self.columns_data = columns_data

    def prepare_column_header(self):
        """Prepare the column header line and finish the header_lines"""
        N_header_lines = len(self.header_lines) + 3
        self.header_lines.append(f"N_header_lines = {N_header_lines}\n")
        self.header_lines.append("\n")

        col_header_line = self.delim.join(self.s_list) + "\n"
        self.header_lines.append(col_header_line)

    def write_header(self):
        """Create the file and write the header lines."""
        with open(self.path_to_file, "w") as f:
            f.writelines(self.header_lines)

    def write_data(self):
        """Write data to the file one line at a time."""
        max_length = max([len(data) for data in self.columns_data.values()])
        for n in range(max_length):
            data_strings = []
            for s_name in self.s_list:
                if len(self.columns_data[s_name]) > n:
                    # Then there's more data to write for this series
                    data_strings.append(str(self.columns_data[s_name][n]))
                else:
                    # Then all this series is written. Just leave space.
                    data_strings.append("")
            line = self.delim.join(data_strings) + "\n"
            with open(self.path_to_file, "a") as f:
                f.write(line)
