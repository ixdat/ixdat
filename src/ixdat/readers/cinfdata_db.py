"""Module defining direct DB reader connection to Surfcat's cinfdata system"""
import warnings
from .. import Measurement
from ..data_series import DataSeries, ValueSeries, TimeSeries, Field
from ..techniques.ms import MSSpectrum
from ..spectra import Spectrum, SpectrumSeries
from ..config import plugins

SCALE_TIME_TO_SECONDS = 1e-3


class CinfdataDBReader:
    """A class that connects to cinf_database or read from cache
    https://cinfdata-dababase-client.readthedocs.io/en/latest/index.html

    Attributes:
        setup_name (str): The setup name in the DB
        timestamp (str): Timestamp when the experiment started in YYYY-MM-DD HH:MM:SS
        units (dict): Dictionary of columns names with corresponding units
        tstamp (str): The unix time corresponding to t=0 for the measurement
        technique (str): The name of the technique
        measurement (Measurement): The measurement returned by read() when the database
            is read. self.measurement is None before read() is called.
    """

    def __init__(self):
        """Initialize a Reader for cinf_database. See class docstring."""
        self.name = None
        self.sample_name = None
        self.setup_name = None
        self.timestamp = None
        self.comment = None
        self.grouping_column = None
        self.tstamp = None
        self.data_has_been_fetched = False
        self.metadata = {}
        self.technique = "reactor"
        self.measurement_class = None
        self.measurement = None
        self.cinf_db = None
        self.include_mass_scans = False
        self.spectrum_list = []
        self.verbose = False

        plugins.activate_cinfdata()

    def read(self, path_to_file, name=None, cls=None, units=None, **kwargs):
        """
        Return a xx-Measurement or Spectrum with the data and metadata recorded from
        a setup at SurfCat at given timestamp
        All attributes of this reader can be accessed from the
        measurement as `measurement.reader.attribute_name`.

        Args:
        path_to_file (str): Named argument from Measurement Class.
            Can be used as the setup name in the cinfdatabase
        **kwargs (dict): Key-word arguments are passed to cinf Measurement.__init__
            setup_name (str): The setup name in the database default to path_to_file
            timestamp (str): Timestamp the measurement started given as
                (YYYY-MM-DD HH:MM:SS)
        """

        if self.data_has_been_fetched:
            print(
                f"This {self.__class__.__name__} has already fetched data from "
                f" {self.token} grouped by {self.grouping_column}"
                " Returning the measurement resulting from the original read. "
                "Use a new Reader if you want to read another file."
            )

            return self.measurement

        self.measurement_class = kwargs.pop("measurement_class", cls)
        self.setup_name = kwargs.pop("setup_name", path_to_file)
        self.timestamp = kwargs.pop("timestamp", None)
        self.comment = kwargs.pop("comment", None)
        self.include_mass_scans = kwargs.pop("include_mass_scans", False)
        self.verbose = kwargs.pop("verbose", False)
        self.grouping_column = kwargs.pop("group", None)

        # figure out whether to collect data group by 'comment' or 'timestamp'

        if not self.grouping_column:
            if self.timestamp and not self.comment:
                self.grouping_column = "time"
                self.token = self.timestamp
            elif self.comment and not self.timestamp:
                self.grouping_column = "comment"
                self.token = self.comment
            else:
                warnings.warn(
                    "Both a comment and a timestamp is given "
                    "but no explicit grouping column is set. \n"
                    f"Defaults to 'time' using '{self.timestamp}'",
                    stacklevel=2,
                )
                self.grouping_column = "time"
                self.token = self.timestamp

        self.data_has_been_fetched = True

        if issubclass(self.measurement_class, Spectrum):
            return self.read_spectrums()

        if issubclass(self.measurement_class, Measurement):
            obj_as_dict = self.read_ms()
            if issubclass(cls, self.measurement_class):
                self.measurement_class = cls
            obj_as_dict.update(kwargs)
            self.measurement = self.measurement_class.from_dict(obj_as_dict)

            if self.include_mass_scans:
                if self.verbose:
                    print("adding mass scans to the measurement")
                self.add_mass_scans()

            return self.measurement

    def read_ms(self):
        """Download MS data from cinfdata_database and make corresponding tseries and
        vseries to place in a dictionary to return.
        return obj_as_dict (dict)
        """

        with plugins.cinfdata(
            setup_name=self.setup_name, grouping_column=self.grouping_column
        ) as cinf_db:

            self.group_data = cinf_db.get_data_group(
                self.token, scaling_factors=(SCALE_TIME_TO_SECONDS, None)
            )

            self.group_meta = cinf_db.get_metadata_group(self.token)

        self.set_sample_name_tstamp_and_name()

        if self.verbose:
            print("Retriving data from measurement named: ", self.sample_name)
            print("Measurement started recording on: ", self.timestamp)

        data_series_list = []

        if self.verbose:
            print("Column names in measurement: ")

        for key in self.group_data.keys():
            meta = self.group_meta[key]
            if meta["type"] != 5:  # 5 is specific mass_time measurements
                continue

            column_name = meta["mass_label"]
            if self.verbose:
                print("Col name: ", column_name)

            tcol = self.group_data[key][:, 0]
            vcol = self.group_data[key][:, 1]

            tseries = TimeSeries(
                name=column_name + "-x",
                unit_name=get_column_unit(column_name + "-x") or "s",
                data=tcol,
                tstamp=self.tstamp,
            )

            vseries = ValueSeries(
                name=column_name,
                data=vcol,
                tseries=tseries,
                unit_name=get_column_unit(column_name + "-y"),
            )
            data_series_list.append(tseries)
            data_series_list.append(vseries)

        obj_as_dict = dict(
            name=self.name,
            sample_name=self.sample_name,
            technique=self.technique,
            reader=self,
            series_list=data_series_list,
            tstamp=self.tstamp,
            metadata=self.group_meta,
        )

        if not data_series_list:
            warnings.warn(
                f"No mass spec data was found using '{self.token}' "
                f" and group_column: '{self.grouping_column}'",
                stacklevel=2,
            )
            return None

        return obj_as_dict

    def read_spectrums(self, **kwargs):
        """
        Download spectrums from cinfdata_database with either a timestamp or a comment.
        When all associated spectras have been downloaded and added to
        self.spectrum_list, this method will try to:
            1. figure out if any spectras were able to be downloaded; if none,return None
            2. figure out if this method was used to add mass_scans to another
               measurement; return list of Spectrum to be added to that measurement
            3. figure out if one and only one spectrum was downloaded; return that single
               Spectrum as the first spectrum in self.spectrum_list[0]
            4. figure out if all spectras downloaded are associated enough to be a
               SpectrumSeries (They have to be equal in dimensions like mass scans)
               (multiple mass scans over time versus XPS spectra from different regions);
               return SpectrumSeries
            5. If not one and only one Spectrum was read and the multiple of spectras
               read were not possible to combine in a SpectrumSeries,
               return self.spectrum_list
        """

        with plugins.cinfdata(
            setup_name=self.setup_name, grouping_column=self.grouping_column
        ) as cinf_db:

            # return dict with measurements as key containing x,y values in a np.array
            self.group_data = cinf_db.get_data_group(self.token)

            # return dict of meta data associated with the key associated (measurement)
            self.group_meta = cinf_db.get_metadata_group(self.token)

        # set sample_name, tstamp and measurement name from meta data from cinfdatabase
        self.set_sample_name_tstamp_and_name()

        self.spectrum_list = []
        for key in self.group_meta:  # key is unique to each measurement
            group_type = self.group_meta[key]["type"]
            if group_type == 2:  # type 2 is unique describing XPS data
                obj_as_dict = self.get_spectrum_as_dict(key, group_type)
                obj_as_dict["name"] = self.group_meta[key]["name"]
                self.spectrum_list.append(self.measurement_class.from_dict(obj_as_dict))
            elif group_type == 4:  # type 4 is unique describing mass scan
                obj_as_dict = self.get_spectrum_as_dict(key, group_type)
                obj_as_dict["name"] = self.sample_name
                self.spectrum_list.append(self.measurement_class.from_dict(obj_as_dict))

        if not self.spectrum_list:
            warnings.warn(
                f"No spectrum was found using '{self.token}' and group_column:"
                f"'{self.grouping_column}'",
                stacklevel=2,
            )
            return None
        elif self.include_mass_scans:
            return self.spectrum_list
        elif len(self.spectrum_list) == 1:
            return self.spectrum_list[0]
        else:
            try:
                return SpectrumSeries.from_spectrum_list(self.spectrum_list)
            except ValueError:
                warnings.warn(
                    "Could not return SpectrumSeries from list of spectrums "
                    f"using '{self.token}' and group column: "
                    f"'{self.grouping_column}'."
                    "\nReturn list of all Spectrums.",
                    stacklevel=2,
                )
                return self.spectrum_list

    def get_spectrum_as_dict(self, key, group_type):
        """Convert spectrum data from cinfdatabase to ixdat dict format with xseries and
        fields.
        Return spectrum_as_dict (dict of xseries and Fields)
        """
        # Extract x and y data columns and timestamp from group data and metadata
        x_col = self.group_data[key][:, 0]
        y_col = self.group_data[key][:, 1]
        metadata = self.group_meta[key]
        tstamp = self.group_meta[key]["unixtime"]

        # Create x DataSeries object with appropriate metadata
        x_series = DataSeries(
            data=x_col,
            name=SPECTRUM_METADATA[group_type]["x_name"],
            unit_name=SPECTRUM_METADATA[group_type]["x_unit_name"],
        )

        # Create y Field obj with appropriate metadata and x DataSeries as its only axis
        y_field = Field(
            data=y_col,
            name=SPECTRUM_METADATA[group_type]["field_name"],
            unit_name=SPECTRUM_METADATA[group_type]["field_unit"],
            axes_series=[x_series],
        )

        # Create dictionary with spectrum object metadata and x-y Field object
        spectrum_as_dict = {
            "sample_name": self.sample_name,
            "technique": SPECTRUM_METADATA[group_type]["technique"],
            "field": y_field,
            "tstamp": tstamp,
            "reader": self,
            "metadata": metadata,
        }

        return spectrum_as_dict

    def add_mass_scans(self):
        """Add corresponding mass scans to mass_time from 'comment'"""
        self.include_mass_scans = True
        self.measurement_class = MSSpectrum
        self.grouping_column = "comment"
        self.token = self.sample_name
        self.spectrum_list = self.read_spectrums()

        if self.verbose:
            print(f"Using {self.measurement.time_series[-1]} to find end of experiment")
            print(
                "Unixtime end of exp "
                f"{self.measurement.time_series[-1].data[-1] + self.tstamp}"
            )

        # Find index of the first spectrum to include as mass scans of this measurement
        first_index = 0
        if self.spectrum_list[0].tstamp < self.tstamp:
            for i, spectrum in enumerate(self.spectrum_list):
                if spectrum.tstamp < self.tstamp:
                    first_index = i
                else:
                    break
        # Find index of the last spectrum to include as mass scans
        last_index = -1
        measurement_end_time = self.measurement.time_series[-1].data[-1] + self.tstamp

        if self.spectrum_list[-1].tstamp > measurement_end_time:
            for i, spectrum in reversed(list(enumerate(self.spectrum_list))):
                if spectrum.tstamp >= measurement_end_time:
                    last_index = i
                else:
                    break

        if self.verbose:
            print(f"Start and End index of spectrum list: {first_index}, {last_index}\n")
            print(
                "Tstamp of first and last spectrum in list: "
                f"{self.spectrum_list[first_index].tstamp}, "
                f"{self.spectrum_list[last_index].tstamp}"
            )

        # Create spectrum series and add it to the measurement
        spectrums_to_add = self.spectrum_list[first_index:last_index]
        ms_spectra = SpectrumSeries.from_spectrum_list(spectrums_to_add)
        self.measurement = self.measurement + ms_spectra

    def set_sample_name_tstamp_and_name(self):
        """Set the sample name and measurement name and tstamp from meta data retrieved
        from cinfdata"""

        metadata = list(self.group_meta.items())[0][1]

        self.sample_name = None
        for key_name in ("Comment", "comment"):
            if key_name in metadata:
                self.sample_name = metadata[key_name]

        if self.sample_name is None:
            print("No comment to set as sample_name.")

        self.name = metadata["time"].strftime("%Y-%m-%d %H:%M:%S")
        self.tstamp = float(metadata["unixtime"])


def get_column_unit(column_name):
    """Return the unit name of an ixdat column, i.e the part of the name after the '/'"""
    if column_name.startswith("M") and column_name.endswith("-y"):
        unit_name = "A"
    elif column_name.startswith("M") and column_name.endswith("-x"):
        unit_name = "s"
    elif column_name.startswith("Reactor") and column_name.endswith("pressure-y"):
        unit_name = "bar"
    elif not column_name.startswith("Reactor") and column_name.endswith("pressure-y"):
        unit_name = "mbar"
    elif column_name.endswith("temperature-y"):
        unit_name = "celcius"
    elif column_name.startswith("Flow"):
        unit_name = "ml/min"

    else:
        # TODO: Figure out how cinfdata represents units for other stuff.
        #    see https://github.com/ixdat/ixdat/pull/30/files#r811432543, and
        #    https://github.com/CINF/cinfdata/blob/master/sym-files2/export_data.py#L125
        unit_name = None

    return unit_name


SPECTRUM_METADATA = {
    2: {
        "x_name": "Binding energy / eV",
        "x_unit_name": "eV",
        "field_name": "Counts per second",
        "field_unit": "n/s",
        "technique": "XPS",
    },
    4: {  # type 4 is unique describing msscan
        "x_name": "Mass [AMU]",
        "x_unit_name": "m/z",
        "field_name": "Current",
        "field_unit": "[A]",
        "technique": "MS_spectrum",
    },
}
