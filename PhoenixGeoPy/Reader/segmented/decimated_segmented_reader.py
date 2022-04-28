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

from struct import unpack_from
from PhoenixGeoPy.Reader import TSReaderBase

# =============================================================================


class DecimatedSegmentedReader(TSReaderBase):
    """Class to create a streamer for segmented decimated time series,
    i.e. *.td_24k"""

    def __init__(self, path, num_files=1, report_hw_sat=False):
        # Init the base class
        super().__init__(path, num_files, 128, report_hw_sat)
        self.unpack_header()
        self.subheader = {}

    def read_subheader(self):
        subheaderBytes = self.stream.read(32)
        if not subheaderBytes:
            if self.open_next():
                subheaderBytes = self.stream.read(32)

        if not subheaderBytes or len(subheaderBytes) < 32:
            self.subheader["timestamp"] = 0
            self.subheader["samplesInRecord"] = 0
        else:
            self.subheader["timestamp"] = unpack_from("I", subheaderBytes, 0)[0]
            self.subheader["samplesInRecord"] = unpack_from("I", subheaderBytes, 4)[0]
            self.subheader["satCount"] = unpack_from("H", subheaderBytes, 8)[0]
            self.subheader["missCount"] = unpack_from("H", subheaderBytes, 10)[0]
            self.subheader["minVal"] = unpack_from("f", subheaderBytes, 12)[0]
            self.subheader["maxVal"] = unpack_from("f", subheaderBytes, 16)[0]
            self.subheader["avgVal"] = unpack_from("f", subheaderBytes, 20)[0]

    def read_record_data(self):
        ret_array = np.empty([0])
        if (
            self.stream is not None
            and self.subheader["samplesInRecord"] is not None
            and self.subheader["samplesInRecord"] != 0
        ):
            ret_array = np.fromfile(
                self.stream, dtype=np.float32, count=self.subheader["samplesInRecord"]
            )
            if ret_array.size == 0:
                if not self.open_next():
                    return np.empty([0])
                # Array below will contain the data, or will be an np.empty array if end of series as desired
                ret_array = np.fromfile(
                    self.stream,
                    dtype=np.float32,
                    count=self.subheader["samplesInRecord"],
                )

        return ret_array

    def read_record(self):
        self.read_subheader()
        return self.read_record_data()
