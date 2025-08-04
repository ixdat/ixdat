# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 11:03:18 2025

@author: Søren
@contributor: Frederik
"""
from ixdat.measurement_base import Measurement
from ixdat.techniques.cv import CyclicVoltammogram


class EChemDBReader:

    def read(self, echemdb_identifier, cls):
        """

        Args:
            echemdb_identifier: the name on https://www.echemdb.org/
            cls: class to return. If cls is Measurement, return a
                CyclicVoltammogram object

        Returns an object of type cls
        """

        if cls is Measurement:
            cls = CyclicVoltammogram
