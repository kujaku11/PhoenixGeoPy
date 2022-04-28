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
    """
    Class to create a streamer for continuous decimated time series,
    i.e. *.td_150, *.td_30
    
    These files have no sub header information.
    """

    def __init__(self, path, num_files=1, report_hw_sat=False):
        # Init the base class
        super().__init__(
            path, num_files=num_files, header_size=128, report_hw_sat=report_hw_sat
            ) 

        self.unpack_header(self.stream)
        self.subheader = {}
    
    # need a read and read sequence
    def read(self):
        """
        Read in the full data from the file given
        
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.stream.seek(self.header_size)
        return np.fromfile(self.stream, dtype=np.float32)
    
    def read_sequence(self, start=0, end=None):
        """
        Read a sequence of files
        
        :param start: DESCRIPTION, defaults to 0
        :type start: TYPE, optional
        :param end: DESCRIPTION, defaults to None
        :type end: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        
        data = np.array([], dtype=np.float32)
        for fn in self.sequence_list[slice(start, end)]:
            self._open_file(fn)
            self.unpack_header(self.stream)
            ts= self.read()
            data = np.append(data, ts)
            
        return data
    
    
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
