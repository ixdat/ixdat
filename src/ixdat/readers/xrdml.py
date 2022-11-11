"""This module defines the reader of .xrdml files from, for example, Empyrion XRD"""

from ixdat import Spectrum
from ixdat.data_series import DataSeries, Field
from xml.dom.minidom import parse
import numpy as np


class XRDMLReader:
    def read(self, path_to_file, cls=None, **kwargs):
        """Read an XRDML file.

        TODO: Finish tutorial here and improve use of xml parser:
          https://realpython.com/python-xml-parser/

        Args:
            path_to_file (str or Path): The path to the .xrdml file to read
            cls (Spectrum class): The class to return an object of.
                Defaults to `Spectrum`.
            kwargs: Additional keyword arguments are passed on to `cls.from_field`
        """
        cls = cls or Spectrum
        with open(path_to_file, "r") as f:
            document = parse(f)
        datapoint_node = document.getElementsByTagName("dataPoints")[0]
        position_nodes = [
            node for node in datapoint_node.childNodes if "positions" in str(node)
        ]
        start_position_element = position_nodes[0].getElementsByTagName("startPosition")[
            0
        ]
        end_position_element = position_nodes[0].getElementsByTagName("endPosition")[0]
        x_min = float(start_position_element.childNodes[0].data)
        x_max = float(end_position_element.childNodes[0].data)
        intensity_parent_node = [
            node
            for node in datapoint_node.childNodes
            if "intensities" in str(node) or "counts" in str(node)
        ][0]
        intensity_node = intensity_parent_node.childNodes[0]
        intensity_string = intensity_node.data
        y_vec = np.fromstring(intensity_string, dtype=float, sep=" ")
        x_vec = np.linspace(x_min, x_max, len(y_vec))
        xseries = DataSeries(name="two theta", unit_name="degree", data=x_vec)
        field = Field(
            name="intensity", unit_name="counts", data=y_vec, axes_series=[xseries]
        )
        return cls.from_field(field, **kwargs)
