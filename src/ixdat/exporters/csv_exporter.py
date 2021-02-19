"""Classes for exporting measurement data"""
from pathlib import Path


class CSVExporter:
    """The default exporter, which writes delimited measurement data row-wise to file"""

    def __init__(self, measurement=None, delim=",\t", default_v_list=None):
        """Initiate the exported with a measurement (Measurement) and delimiter (str)"""
        self.measurement = measurement
        self.delim = delim
        self._default_v_list = default_v_list

    @property
    def default_v_list(self):
        """This will typically be overwritten by inheriting Exporters"""
        return self._default_v_list

    def export(self, *args, **kwargs):
        """Export the exporter's measurement via exporter.export_measurement()"""
        return self.export_measurement(self.measurement, *args, **kwargs)

    def export_measurement(self, measurement, path_to_file, v_list=None, tspan=None):
        """Export a given measurement to a specified file.

        Args:
            measurement (Measurement): The measurement to export
            path_to_file (Path): The path to the file to measure. If it has no suffix,
                a .csv suffix is appended.
            v_list (list of str): The names of the data series to include. Defaults in
                CSVExporter to all VSeries and TSeries in the measurement. This default
                may be overwritten in inheriting exporters.
            tspan (timespan): The timespan to include in the file, defaults to all of it
        """
        columns_data = {}
        s_list = []
        v_list = v_list or self.default_v_list or list(measurement.data_cols)
        if isinstance(path_to_file, str):
            path_to_file = Path(path_to_file)
        if not path_to_file.suffix:
            path_to_file = path_to_file.with_suffix(".csv")

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

        N_header_lines = len(header_lines) + 3
        header_lines.append(f"N_header_lines = {N_header_lines}\n")
        header_lines.append("\n")

        col_header_line = (
            "".join([s_name + self.delim for s_name in s_list])[: -len(self.delim)]
            + "\n"
        )
        header_lines.append(col_header_line)

        lines = header_lines
        max_length = max([len(data) for data in columns_data.values()])
        for n in range(max_length):
            line = ""
            for s_name in s_list:
                if len(columns_data[s_name]) > n:
                    line = line + str(columns_data[s_name][n]) + self.delim
                else:
                    line = line + self.delim
            line = line + "\n"
            lines.append(line)

        with open(path_to_file, "w") as f:
            f.writelines(lines)
