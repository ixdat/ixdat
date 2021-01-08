"""Classes for exporting measurement data"""


class CSVExporter:
    """The default exporter, which writes delimited measurement data row-wise to file"""

    def __init__(self, measurement=None, delim=",\t"):
        """Initiate the exported with a measurement (Measurement) and delimiter (str)"""
        self.measurement = measurement
        self.delim = delim

    def export(self, *args, **kwargs):
        """Export the exporter's measurement via exporter.export_measurement()"""
        return self.export_measurement(self.measurement, *args, **kwargs)

    def export_measurement(self, measurement, path_to_file, v_list=None, tspan=None):
        """Export a given measurement to a specified file.

        Args:
            measurement (Measurement): The measurement to export
            path_to_file (Path): The path to the file to measure. If it has no suffix,
                a .csv suffix is appended.
            v_list (list of str): The names of the data series to include. Defaults to
                all VSeries and TSeries in the measurement.
            tspan (timespan): The timespan to include in the file, defaults to all of it
        """
        columns_data = {}
        s_list = []
        v_list = v_list or list(measurement.data_cols)
        if not path_to_file.metadata_suffix:
            path_to_file = path_to_file.with_suffix(".csv")
        for v_name in v_list:
            t_name = measurement[v_name].tseries.name
            t, v = measurement.get_t_and_v(v_name, tspan=tspan)
            if t_name not in columns_data:
                columns_data[t_name] = t
                s_list.append(t_name)
            columns_data[v_name] = v
            s_list.append(v_name)

        header_line = (
            "".join([s_name + self.delim for s_name in s_list])[: -len(self.delim)]
            + "\n"
        )

        lines = [header_line]
        max_length = max([len(data) for data in columns_data.values()])
        for n in range(max_length):
            line = ""
            for s_name in s_list:
                if len(columns_data[s_name]) > n:
                    line.join(str(columns_data[s_name][n]) + self.delim)
                else:
                    line.join(self.delim)
            lines.append(line)

        with open(path_to_file, "w") as f:
            f.writelines(lines)
