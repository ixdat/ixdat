import numpy as np
from .db import Saveable, PlaceHolderObject
from .data_series import DataSeries, TimeSeries, Field


class Spectrum(Saveable):
    """The Spectrum class"""

    table_name = "spectrum"
    column_attrs = {
        "name",
        "technique",
        "metadata",
        "sample_name",
        "tstamp",
        "y_id",
    }

    def __init__(
        self,
        name,
        technique=None,
        metadata=None,
        sample_name=None,
        reader=None,
        field=None,
        field_id=None,
    ):
        """
        Args:
            name (str): The name of the spectrum
            metadata (dict): Free-form spectrum metadata. Must be json-compatible.
            technique (str): The spectrum technique
            sample_name (str): The sample name
            reader (Reader): The reader, if read from file
            field (Field): The Field containing the data (x, y, and tstamp)
            field_id (id): The id in the data_series table of the Field with the data, if
                the field is not yet loaded from backend.
        """
        super().__init__()
        self.name = name
        self.technique = technique
        self.metadata = metadata
        self.sample_name = sample_name
        self.reader = reader
        self._field = field or PlaceHolderObject(field_id, cls=Field)

        self._plotter = None
        # defining this method here gets it the right docstrings :D
        self.plot = self.plotter.plot

    @property
    def plotter(self):
        """The default plotter for Measurement is ValuePlotter."""
        if not self._plotter:
            from .plotters.spectrum_plotter import SpectrumPlotter

            # FIXME: I had to import here to avoid running into circular import issues

            self._plotter = SpectrumPlotter(spectrum=self)
        return self._plotter

    @classmethod
    def read(cls, path_to_file, reader, **kwargs):
        """Return a Measurement object from parsing a file with the specified reader

        Args:
            path_to_file (Path or str): The path to the file to read
            reader (str or Reader class): The (name of the) reader to read the file with.
            kwargs: key-word arguments are passed on to the reader's read() method.
        """
        if isinstance(reader, str):
            # TODO: see if there isn't a way to put the import at the top of the module.
            #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437471
            from .readers import READER_CLASSES

            reader = READER_CLASSES[reader]()
        # print(f"{__name__}. cls={cls}")  # debugging
        return reader.read(path_to_file, cls=cls, **kwargs)

    @classmethod
    def from_data(
        cls,
        x,
        y,
        tstamp=None,
        x_name="x",
        y_name="y",
        x_unit_name=None,
        y_unit_name=None,
        **kwargs
    ):
        xseries = DataSeries(data=x, name=x_name, unit_name=x_unit_name)
        yseries = DataSeries(data=y, name=y_name, unit_name=y_unit_name)
        return cls.from_series(xseries, yseries, tstamp, **kwargs)

    @classmethod
    def from_series(cls, xseries, yseries, tstamp, **kwargs):
        tseries = TimeSeries(
            data=np.array([0]), tstamp=tstamp, unit_name="s", name="spectrum time / [s]"
        )
        field = Field(
            data=yseries.data,
            axes_series=[xseries, tseries],
            name=yseries.name,
            unit_name=yseries.unit_name,
        )
        return cls.from_field(field, **kwargs)

    @classmethod
    def from_field(cls, field, **kwargs):
        spectrum_as_dict = kwargs
        spectrum_as_dict["field"] = field
        if "name" not in spectrum_as_dict:
            spectrum_as_dict["name"] = field.name
        return cls.from_dict(spectrum_as_dict)

    @property
    def field(self):
        if isinstance(self._field, PlaceHolderObject):
            self._field = self._field.get_object()
        return self._field

    @property
    def xseries(self):
        return self.field.axes_series[0]

    @property
    def x(self):
        return self.xseries.data

    @property
    def x_name(self):
        return self.xseries.name

    @property
    def yseries(self):
        return DataSeries(
            name=self.field.name, data=self.field.data, unit_name=self.field.unit_name
        )

    @property
    def y(self):
        return self.field.data

    @property
    def y_name(self):
        return self.field.name

    @property
    def tseries(self):
        return self.field.axes_series[1]

    @property
    def tstamp(self):
        tseries = self.tseries
        return tseries.data[0] + tseries.tstamp
