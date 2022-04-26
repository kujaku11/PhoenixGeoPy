# -*- coding: utf-8 -*-
"""Module to read and parse native Phoenix Geophysics data formats of the MTU-5C Family

This module implements Streamed readers for segmented-decimated continuus-decimated
and native sampling rate time series formats of the MTU-5C family.
"""

__author__ = "Jorge Torres-Solis"

# =============================================================================
# Imports
# =============================================================================

import numpy as np

from PhoenixGeoPy.Reader import TSReaderBase

# =============================================================================


class DecimatedContinuousReader(TSReaderBase):
    """Class to create a streamer for continuous decimated time series,
    i.e. *.td_150, *.td_30"""

    def __init__(self, path, num_files=1, report_hw_sat=False):
        # Init the base class
        super().__init__(path, num_files, 128, report_hw_sat)
        self.unpack_header()
        self.subheader = {}

    def read_data(self, numSamples):
        ret_array = np.empty([0])
        if self.stream is not None:
            ret_array = np.fromfile(self.stream, dtype=np.float32, count=numSamples)
            while ret_array.size < numSamples:
                if not self.open_next():
                    return np.empty([0])
                # Array below will contain the data, or will be an np.empty array if end of series as desired
                ret_array = np.append(
                    ret_array,
                    np.fromfile(
                        self.stream,
                        dtype=np.float32,
                        count=(numSamples - ret_array.size),
                    ),
                )
        return ret_array
