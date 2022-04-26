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

from struct import unpack_from, unpack
from PhoenixGeoPy.Reader import DataScaling, TSReaderBase

# =============================================================================


class NativeReader(TSReaderBase):
    """Native sampling rate 'Raw' time series reader class"""

    def __init__(
        self,
        path,
        num_files=1,
        scale_to=DataScaling.AD_input_volts,
        header_size=128,
        last_frame=0,
        ad_plus_minus_range=5.0,
        channel_type="E",
        report_hw_sat=False,
    ):
        # Init the base class
        super().__init__(path, num_files, header_size, report_hw_sat)

        self._chunk_size = 4096

        # Track the last frame seen by the streamer, to report missing frames
        self.last_frame = last_frame
        self.header_size = header_size
        self.data_scaling = scale_to
        self.ad_plus_minus_range = ad_plus_minus_range

        if header_size == 128:
            self.unpack_header(self.stream)

        # Now that we have the channel circuit-based gain (either form init or from the header)
        # We can calculate the voltage range at the input of the board.
        self.input_plusminus_range = self._calculate_input_plusminus_range()

        self.scale_factor = self._calculate_data_scaling()

        # optimization variables
        self.footer_idx_samp_mask = int("0x0fffffff", 16)
        self.footer_sat_mask = int("0x70000000", 16)

    def _calculate_input_plusminus_range(self):
        """
        set the correct input plusminus range from metadata
        """
        return self.ad_plus_minus_range / self.total_circuitry_gain

    def _calculate_data_scaling(self):
        """
        Get the correct data scaling for the AD converter
        """
        if self.data_scaling == DataScaling.AD_in_ADunits:
            return 256
        elif self.data_scaling == DataScaling.AD_input_volts:
            return self.ad_plus_minus_range / (2 ** 31)
        elif self.data_scaling == DataScaling.instrument_input_volts:
            return self.input_plusminus_range / (2 ** 31)
        else:
            raise LookupError("Invalid scaling requested")

    def read_frames(self, num_frames):
        """
        Read the given amount of frames from the data.

        .. note:: that seek is not reset so if you iterate this the stream
        reads from the last tell.

        :param num_frames: DESCRIPTION
        :type num_frames: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        frames_in_buf = 0
        _idx_buf = 0
        _data_buf = np.empty([num_frames * 20])  # 20 samples packed in a frame

        while frames_in_buf < num_frames:

            dataFrame = self.stream.read(64)
            if not dataFrame:
                if not self.open_next():
                    return np.empty([0])
                dataFrame = self.stream.read(64)
                if not dataFrame:
                    return np.empty([0])

            dataFooter = unpack_from("I", dataFrame, 60)

            # Check that there are no skipped frames
            frameCount = dataFooter[0] & self.footer_idx_samp_mask
            difCount = frameCount - self.last_frame
            if difCount != 1:
                print(
                    "Ch [%s] Missing frames at %d [%d]\n"
                    % (self.channel_id, frameCount, difCount)
                )
            self.last_frame = frameCount

            for ptrSamp in range(0, 60, 3):
                # unpack expects 4 bytes, but the frames are only 3?
                value = unpack(">i", dataFrame[ptrSamp : ptrSamp + 3] + b"\x00")[0]
                _data_buf[_idx_buf] = value * self.scale_factor
                _idx_buf += 1

            frames_in_buf += 1

            if self.report_hw_sat:
                satCount = (dataFooter[0] & self.footer_sat_mask) >> 24
                if satCount:
                    print(
                        "Ch [%s] Frame %d has %d saturations"
                        % (self.ch_id, frameCount, satCount)
                    )

        return _data_buf

    @property
    def npts_per_frame(self):
        return int((self.frame_size_bytes - 4) / 3)

    def _get_number_of_frames(self):
        return int((self.file_size - self.header_size) / self.frame_size_bytes)

    def _get_number_of_chunks(self):
        # will need to take into account residual if the chunk size is not
        # a good choice.
        return int((self.file_size - self.header_size) / self._chunk_size)
    
    def _get_number_of_frames_per_chunk(self):
        return int((self._chunk_size / self.frame_size_bytes) * self.npts_per_frame)

    def read(self):

        # first read in whole file, probably should do this in chunks

        n_frames = self._get_number_of_frames()
        data = np.zeros(n_frames * self.npts_per_frame, dtype=np.int32)
        footer = np.zeros(n_frames)

        # start from the end of the header
        self.stream.seek(self.header_size)

        chunk_frames = int(self._chunk_size / self.frame_size_bytes)
        data_slices = [
            slice(ii * 64, (ii + 1) * 60 + 4 * ii) for ii in range(chunk_frames)
        ]
        footer_slices = [
            slice((ii) * 60, (ii) * 60 + 4) for ii in range(1, chunk_frames + 1)
        ]
        for count in range(self._get_number_of_chunks()):
            byte_string = self.stream.read(self._chunk_size)
            
            # maybe try to read it all into a numpy array and then reshape
            # add the 4 byte and convert?
            # np.frombuffer(s.read(), dtype=np.dtype("u1"))
            # reshape to (nframes, 64)
            # split footer off c = b[:, 0:60].
            # f = int(c.shape / 12)
            # f = int(c.size / 12)
            # f
            # from numpy.lib.stride_tricks import as_strided
            # rd = as_strided(c.view(np.int32), strides=(12,3), shape=(f, 4))

            for data_frame, footer_frame, ii in zip(
                data_slices, footer_slices, range(n_frames)
            ):
                # need to make this part more efficient, should use numpy 
                # for this, should be faster and wouldn't have to loop
                index_0 = (count * chunk_frames + ii) * self.npts_per_frame
                index_1 = (count * chunk_frames + ii + 1) * self.npts_per_frame
                data[index_0:index_1] = [
                    unpack(
                        ">i",
                        byte_string[data_frame][slice(jj * 3, (jj + 1) * 3)] + b"\x00",
                    )[0]
                    for jj in range(self.npts_per_frame)
                ]
                footer[ii] = unpack("I", byte_string[footer_frame])[0]

        return data, footer

    def read_sequence(self, start=0, end=None):
        """
        Read sequences

        :param start: DESCRIPTION, defaults to 0
        :type start: TYPE, optional
        :param end: DESCRIPTION, defaults to None
        :type end: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        data = np.array([])
        for fn in self.sequence_list[slice(start, end)]:
            self._open_file(fn)
            self.unpack_header(self.stream)
            data = np.append(data, self.read())

        return data

    def skip_frames(self, num_frames):
        bytes_to_skip = int(num_frames * 64)
        # Python is dumb for seek and tell, it cannot tell us if a seek goes
        # past EOF so instead we need to do inefficient reads to skip bytes
        while bytes_to_skip > 0:
            foo = self.stream.read(bytes_to_skip)
            local_read_size = len(foo)

            # If we ran out of data in this file before finishing the skip,
            # open the next file and return false if there is no next file
            # to indicate that the skip ran out of
            # data before completion
            if local_read_size == 0:
                more_data = self.open_next()
                if not more_data:
                    return False
            else:
                bytes_to_skip -= local_read_size

        # If we reached here we managed to skip all the data requested
        # return true
        self.last_frame += num_frames
        return True
