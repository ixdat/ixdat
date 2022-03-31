"""Module defining direct DB reader connection to Surfcat's legendary cinfdata system"""

from pathlib import Path
from cinfdata import Cinfdata
import numpy as np
from ..exceptions import ReadError
from ..data_series import ValueSeries, TimeSeries
from ..techniques import MSMeasurement
from .reading_tools import timestamp_string_to_tstamp


class CinfDatabaseReader:
    """A class that connects to cinf_database or read from cache
    https://cinfdata-dababase-client.readthedocs.io/en/latest/index.html

    Attributes:
        setup_name (str): The setup name in the DB
        timestamp (str): Timestamp when the experiment started in YYYY-MM-DD HH:MM:SS
        units (dict): Dictionary of columns names with corrosponding units
        tstamp (str): The unix time corresponding to t=0 for the measurement

        tstamp_list (list of float): list of epoch tstamps in the file's timestamp line

        column_tstamps (dict): The unix time corresponding to t=0 for each time column
        technique (str): The name of the technique
        column_names (list of str): The names of the data columns in the database
        t_and_v_cols (dict): {name: (tcol, vcol)} where name is the name of the
            ValueSeries (e.g. "M2"), tcol is the name of the corresponding time column
            in the file (e.g. "M2-x"), and vcol is the the name of the value column in
            the file (e.g. "M2-y).
        column_data (dict of str: np.array): The data in the file as a dict.
            Note that the np arrays are the same ones as in the measurement's DataSeries,
            so this does not waste memory.
        measurement (Measurement): The measurement returned by read() when the database is
            read. self.measurement is None before read() is called.
    """

    def __init__(self):
        """Initialize a Reader for cinf_database. See class docstring."""
        self.name = None
        self.sample_name = None
        self.setup_name = None
        self.timestamp = None
        self.tstamp = None
        self.tstamp_list = []
        self.column_tstamps = {}
        self.column_names = []
        self.t_and_v_cols = {}
        self.column_data = {}
        self.data_has_been_fetch = False
        self.metadata = {}
        self.technique = "MS"  # TODO: MS? Figure out how to tell if it's something else
        self.measurement_class = MSMeasurement
        self.measurement = None
        self.db = None

    def read(self, setup_name, timestamp, name=None, cls=None, units=None, **kwargs):
        """Return a MSMeasurement with the data and metadata recorded from 
        a setup at SurfCat at given timestamp

        MSMeasurement contains a reference to the reader.

        All attributes of this reader can be accessed from the
        measurement as `measurement.reader.attribute_name`.

        Args:
            setup (str): The setup name in the database
            timestamp (str): Timestamp the measurement started (YYYY-MM-DD HH:MM:SS)
            **kwargs (dict): Key-word arguments are passed to cinf Measurement.__init__
        """

        self.db = Cinfdata(setup_name=setup_name, grouping_column = 'time')

        self.group_data = self.db.get_data_group(timestamp, scaling_factors=(1E-3, None))
        self.group_meta = self.db.get_metadata_group(timestamp, scaling_factors=(1E-3, None))
        self.meta = self.group_meta[list(self.group_meta.keys())[0]]  #



        for key in self.group_data.keys():
            column_name = self.group_meta[key]['mass_label']
            tstamp = self.group_meta[key]['unixtime']
            tcol  = self.group_data[key][:, 0]
            vcol = self.group_data[key][:, 1]

            tseries = TimeSeries(
                    name=column_name+"-x",
                    unit_name=get_column_unit(name) or "s",
                    data=tcol,
                    tstamp=tstamp,
            )
            vseries = ValueSeries(
                    name=column_name,
                    data=vcol,
                    tseries=tseries,
                    unit_name=get_column_unit(vcol),
            )
            data_series_list.append(tseries)
            data_series_list.append(vseries)

        obj_as_dict = dict(
            name=self.name,
            technique=self.technique,
            reader=self,
            series_list=data_series_list,
            tstamp=self.tstamp,
        )


        obj_as_dict.update(kwargs)

        if issubclass(cls, self.measurement_class):
            self.measurement_class = cls

        self.measurement = self.measurement_class.from_dict(obj_as_dict)
        self.data_has_been_fetch = True
        return self.measurement

    def set_sample_name(self):
        try:
            self.sample_name = self.meta['Comment']
        except KeyError:
            try:
                self.sample_name = self.meta['comment']
            except KeyError as e:
                self.sample_name = None
                print('No comment to set as sample_name. ', e)

    def set_name(self):
        self.name = self.meta['time'].strftime("%Y-%m-%d %H:%M:%S")

    def set_tstamp(self):
        self.tstamp = meta['unixtime']


    def add_mass_scans(self):
        ''' Get corrosponding mass scans to mass_time from 'comment' '''
        db = Cinfdata(setup=self.setup, grouping_column='comment')
        group_meta = db.get_metadata_group(self.sample_name)
        group_data = db.get_data_group(self.sample_name)

        for i, key in enumerate(group_meta.keys()):
            if group_meta[key]['mass_label'] == 'Mass Scan' or \
            group_meta[key]['type'] == 4:
                pass  #create spectrum with unixtime to stich in in the measruement
            else:
                pass


def get_column_unit(column_name):
    """Return the unit name of an ixdat column, i.e the part of the name after the '/'"""
    if column_name.startswith("M") and column_name.endswith("-y"):
        unit_name = "A"
    elif column_name.startswith("M") and column_name.endswith("-x"):
        unit_name = "s"
    else:
        # TODO: Figure out how cinfdata represents units for other stuff.
        #    see https://github.com/ixdat/ixdat/pull/30/files#r811432543, and
        #    https://github.com/CINF/cinfdata/blob/master/sym-files2/export_data.py#L125
        unit_name = None
    return unit_name
