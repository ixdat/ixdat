"""Classes for exporting measurement data"""
from pathlib import Path


class CSVExporter:
    """The default exporter, which writes delimited measurement data row-wise to file"""

    default_v_list = None  # This will typically be overwritten by inheriting Exporters

    def __init__(self, measurement=None, delim=",\t"):
        """Initiate the exported with a measurement (Measurement) and delimiter (str)"""
        self.measurement = measurement
        self.delim = delim
        self.header_lines = None
        self.s_list = None
        self.columns_data = None
        self.path_to_file = None

    def export(self, measurement=None, path_to_file=None, v_list=None, tspan=None):
        """Export a given measurement to a specified file.

        To improve flexibility with inheritance, this method allocates its work to:
        - CSVExporter.prepare_header_and_data()
        - CSVExporter.write_header()
        - CSVExporter.write_data()

        Args:
            measurement (Measurement): The measurement to export.
                Defaults to self.measurement.
            path_to_file (Path): The path to the file to write. If it has no suffix,
                a .csv suffix is appended. Defaults to "{measurement.name}.csv"
            v_list (list of str): The names of the data series to include. Defaults in
                CSVExporter to all VSeries and TSeries in the measurement. This default
                may be overwritten in inheriting exporters.
            tspan (timespan): The timespan to include in the file, defaults to all of it
        """
        measurement = measurement or self.measurement
        if not path_to_file:
            path_to_file = input("enter name of file to export.")
        if isinstance(path_to_file, str):
            path_to_file = Path(path_to_file)
        if not path_to_file.suffix:
            path_to_file = path_to_file.with_suffix(".csv")
        self.path_to_file = path_to_file
        self.prepare_header_and_data(measurement, v_list, tspan)
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
        s_list = []
        v_list = v_list or self.default_v_list or list(measurement.value_names)

        timecols = {}
        for v_name in v_list:
            t_name = measurement[v_name].tseries.name
            t, v = measurement.grab(v_name, tspan=tspan)
            if t_name in timecols:
                timecols[t_name].append(v_name)
            else:
                columns_data[t_name] = t
                s_list.append(t_name)
                timecols[t_name] = [v_name]
            columns_data[v_name] = v
            s_list.append(v_name)

        header_lines = []
        for attr in ["name", "technique", "tstamp", "backend_name", "id"]:
            line = f"{attr} = {getattr(measurement, attr)}\n"
            header_lines.append(line)
        for t_name, v_names in timecols.items():
            line = (
                f"timecol '{t_name}' for: "
                + " and ".join([f"'{v_name}'" for v_name in v_names])
                + "\n"
            )
            header_lines.append(line)
        self.header_lines = header_lines
        self.s_list = s_list
        self.columns_data = columns_data

    def prepare_column_header(self):
        """Prepare the column header line and finish the header_lines"""
        N_header_lines = len(self.header_lines) + 3
        self.header_lines.append(f"N_header_lines = {N_header_lines}\n")
        self.header_lines.append("\n")

        col_header_line = (
            "".join([s_name + self.delim for s_name in self.s_list])[: -len(self.delim)]
            + "\n"
        )
        self.header_lines.append(col_header_line)

    def write_header(self):
        """Create the file and write the header lines."""
        with open(self.path_to_file, "w") as f:
            f.writelines(self.header_lines)

    def write_data(self):
        """Write data to the file one line at a time."""
        max_length = max([len(data) for data in self.columns_data.values()])
        for n in range(max_length):
            line = ""
            for s_name in self.s_list:
                if len(self.columns_data[s_name]) > n:
                    # Then there's more data to write for this series
                    line = line + str(self.columns_data[s_name][n]) + self.delim
                else:
                    # Then all this series is written. Just leave space.
                    line = line + self.delim
            line = line + "\n"
            with open(self.path_to_file, "a") as f:
                f.write(line)
