class CSVExporter:
    def __init__(self, measurement=None, delim=",\t"):
        self.measurement = measurement
        self.delim = delim

    def export(self, *args, **kwargs):
        return self.export_measurement(self.measurement, *args, **kwargs)

    def export_measurement(self, measurement, path_to_file, v_list, tspan=None):
        columns_data = {}
        s_list = []
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
