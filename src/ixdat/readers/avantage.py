from pathlib import Path
import numpy as np
from ..spectra import Spectrum
from ..data_series import DataSeries, Field


class AvantageAVGReader:
    """A class for importing a .avg file exported by DataSpace_BatchDump.exe"""

    def __init__(self, path_to_file=None):
        self.path_to_file = path_to_file

    def read(self, path_to_file, cls=None, **kwargs):
        """Load data stored as text by Advantage's default exporting mode

        Copied from pyThetaProbe, written by Anna Winiwarter and Soren Scott in 2019
        TODO: Improve this code. See suggestions here:
                https://github.com/ixdat/ixdat/pull/73#discussion_r892233369
            See also more powerful reader ideas here:
                https://github.com/CINF/PyExpLabSys/blob/master/PyExpLabSys/file_parsers/avantage.py#L384
            and here:
                https://github.com/ixdat/LowOverpotentialRegime/blob/main/src/pyOER/iss.py
        Written for simple intensity-vs-energy, but with possible future expansion in
        mind.
        Returns the dataset as a python dictionary.

        Args:
            path_to_file (str or Path): Path to the .avg data
            cls (Spectrum subclass): Class of spectrum to return an object of
            kwargs: Additional keyword arguments are passed to cls.__init__
        """

        path_to_file = Path(path_to_file or self.path_to_file)
        cls = cls or Spectrum

        get_data_ax_keys, get_data_ax_vals, get_space_ax_keys, get_space_ax_vals = (
            False,
            False,
            False,
            False,
        )
        header_lines = []

        with open(path_to_file, "r") as f:
            in_header = True
            while in_header:
                line = f.readline()
                header_lines += [line]

                if get_data_ax_keys:
                    data_ax_key_str = line.split("=")[-1]
                    data_ax_keys = [k.strip() for k in data_ax_key_str.split(",")]
                    data_axes = {}
                    get_data_ax_keys = False

                if get_space_ax_keys:
                    space_ax_key_str = line.split("=")[-1]
                    space_ax_keys = [k.strip() for k in space_ax_key_str.split(",")]
                    space_axes = {}
                    get_space_ax_keys = False

                if get_data_ax_vals:
                    try:
                        data_ax_nr = int(line.split("=")[0])
                    except ValueError:
                        get_data_ax_vals = False
                        continue
                    data_ax_val_str = line.split("=")[-1]
                    data_ax_vals = [k.strip() for k in data_ax_val_str.split(",")]
                    data_axes[data_ax_nr] = dict(zip(data_ax_keys, data_ax_vals))

                if get_space_ax_vals:
                    try:
                        space_ax_nr = int(line.split("=")[0])
                    except ValueError:
                        get_space_ax_vals = False
                        continue
                    space_ax_val_str = line.split("=")[-1]
                    space_ax_vals = [k.strip() for k in space_ax_val_str.split(",")]
                    space_axes[space_ax_nr] = dict(zip(space_ax_keys, space_ax_vals))

                if "$DATA=" in line:
                    in_header = False
                    in_data = True
                if "data ax" in line:
                    get_data_ax_keys = True
                if "$DATAAXES" in line:
                    get_data_ax_vals = True
                if "space ax" in line:
                    get_space_ax_keys = True
                if "$SPACEAXES" in line:
                    get_space_ax_vals = True

            y_vec = np.array([])
            while in_data:
                line = f.readline()
                data_str = line.split("=")[-1]
                data_str = data_str.replace("#empty#", "nan")
                try:
                    y_i = np.array([float(y_str) for y_str in data_str.split(",")])
                except ValueError:
                    if len(line.strip()) > 0:
                        print(
                            "found no data on this line: \n"
                            + line
                            + "\n Ending data import!"
                        )
                    in_data = False
                else:
                    y_vec = np.append(y_vec, y_i)

            # looks from start and finish values like, for XPS, data axis refers to
            # kinetic energy and space axis refers to binding energy

            x_start = float(space_axes[0]["start"])
            x_width = float(space_axes[0]["width"])
            x_numPoints = int(space_axes[0]["numPoints"])

            x_vec = x_start + np.arange(0, x_numPoints) * x_width

        xseries = DataSeries(
            name=space_axes[0]["label"], unit_name=space_axes[0]["unit"], data=x_vec
        )
        field = Field(name="counts", unit_name="", data=y_vec, axes_series=[xseries])

        if "name" not in kwargs:
            kwargs["name"] = path_to_file.stem

        return cls.from_field(field, **kwargs)
