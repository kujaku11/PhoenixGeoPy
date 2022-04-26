"""Module to read and parse native Phoenix Geophysics data formats of the MTU-5C Family

This module implements Streamed readers for segmented-decimated continuus-decimated
and native sampling rate time series formats of the MTU-5C family.
"""

__author__ = "Jorge Torres-Solis"

# =============================================================================
# Imports
# =============================================================================

from pathlib import Path
import numpy as np

from struct import unpack_from, unpack
from PhoenixGeoPy.Reader import DataScaling, Header

# =============================================================================


class _TSReaderBase(Header):
    def __init__(self, path, num_files=1, header_size=128, report_hw_sat=False):
        self._seq = None
        super().__init__(**{"header_size": header_size, "report_hw_sat": report_hw_sat})
        self.base_path = path
        self.last_seq = self.seq + num_files
        self.stream = None
        # Open the file passed as the fisrt file in the sequence to stream
        self._open_file(self.base_path)

    @property
    def base_path(self):
        return self._base_path

    @base_path.setter
    def base_path(self, value):
        self._base_path = Path(value)

    @property
    def base_dir(self):
        return self.base_path.parent

    @property
    def file_name(self):
        return self.base_path.name

    @property
    def file_extension(self):
        return self.base_path.suffix

    @property
    def inst_id(self):
        return self.base_path.stem.split("_")[0]

    @property
    def rec_id(self):
        if self._rec_id is None:
            return self.base_path.stem.split("_")[1]
        else:
            return self._rec_id

    @rec_id.setter
    def rec_id(self, value):
        self._rec_id = value

    @property
    def ch_id(self):
        if self._ch_id is None:
            return self.base_path.stem.split("_")[2]
        else:
            return self._ch_id

    @ch_id.setter
    def ch_id(self, value):
        self._ch_id = value

    @property
    def seq(self):
        if self._seq is None:
            return int(self.base_path.stem.split("_")[3], 16)
        return self._seq

    @seq.setter
    def seq(self, value):
        self._seq = int(value)

    @property
    def file_size(self):
        return self.base_path.stat().st_size

    @property
    def max_samples(self):
        return int((self.file_size - self.header_size * 3) / 3)

    @property
    def sequence_list(self):
        """
        get all the files in the sequence
        """
        return sorted(list(self.base_dir.glob(f"*{self.file_extension}")))

    def _open_file(self, filename):
        if filename.exists():
            print(f"INFO: Opening {filename}")
            self.stream = open(filename, "rb")
            return True
        return False

    def open_next(self):
        if self.stream is not None:
            self.stream.close()
        self.seq += 1
        self.open_file_seq(self.seq)
        if self.seq < self.last_seq:
            new_path = self.sequence_list[self.seq - 1]
            return self._open_file(new_path)
        return False

    def open_file_seq(self, file_seq_num=None):
        if self.stream is not None:
            self.stream.close()
        if file_seq_num is not None:
            self.seq = file_seq_num
        new_path = self.sequence_list[self.seq - 1]
        return self._open_file(new_path)

    def close(self):
        if self.stream is not None:
            self.stream.close()


class NativeReader(_TSReaderBase):
    """Native sampling rate 'Raw' time series reader class"""

    def __init__(
        self,
        path,
        num_files=1,
        scale_to=DataScaling.AD_input_volts,
        header_size=128,
        last_frame=0,
        channel_gain=0.5,
        ad_plus_minus_range=5.0,
        channel_type="E",
        report_hw_sat=False,
    ):
        # Init the base class
        super().__init__(path, num_files, header_size, report_hw_sat)

        # Track the last frame seen by the streamer, to report missing frames
        self.last_frame = last_frame
        self.header_size = header_size
        self.data_scaling = scale_to
        self.total_circuitry_gain = channel_gain
        self.ad_plus_minus_range = ad_plus_minus_range

        if header_size == 128:
            self.unpack_header(self.stream)

        # Now that we have the channel circuit-based gain (either form init or from the header)
        # We can calculate the voltage range at the input of the board.
        self.input_plusminus_range = (
            self.ad_plus_minus_range / self.total_circuitry_gain
        )

        if self.data_scaling == DataScaling.AD_in_ADunits:
            self._scale_factor = 256
        elif self.data_scaling == DataScaling.AD_input_volts:
            self._scale_factor = self.ad_plus_minus_range / (2 ** 31)
        elif self.data_scaling == DataScaling.instrument_input_volts:
            self._scale_factor = self.input_plusminus_range / (2 ** 31)
        else:
            raise LookupError("Invalid scaling requested")

        # optimization variables
        self.footer_idx_samp_mask = int("0x0fffffff", 16)
        self.footer_sat_mask = int("0x70000000", 16)

    def read_frames(self, num_frames):
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
                    % (self.ch_id, frameCount, difCount)
                )
            self.last_frame = frameCount

            for ptrSamp in range(0, 60, 3):
                # unpack expectes 4 bytes, but the frames are only 3?
                value = unpack(">i", dataFrame[ptrSamp : ptrSamp + 3] + b"\x00")[0]
                _data_buf[_idx_buf] = value * self._scale_factor
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

    def read(self):

        # first read in whole file
        self.stream.seek(self.header_size)
        byte_string = self.stream.read()

        n_frames = int(len(byte_string) / self.frame_size_bytes)
        data_slices = [slice(ii * 64, (ii + 1) * 60 + 4 * ii) for ii in range(n_frames)]
        footer_slices = [
            slice((ii) * 60, (ii) * 60 + 4) for ii in range(1, n_frames + 1)
        ]

        data = np.zeros(n_frames * self.npts_per_frame, dtype=np.int32)
        footer = np.zeros(n_frames)
        for data_frame, footer_frame, ii in zip(
            data_slices, footer_slices, range(n_frames)
        ):
            # need to make this part more efficient
            data[ii * self.npts_per_frame : (ii + 1) * self.npts_per_frame] = [
                unpack(
                    ">i", byte_string[data_frame][slice(jj * 3, (jj + 1) * 3)] + b"\x00"
                )[0]
                for jj in range(self.npts_per_frame)
            ]
            footer[ii] = unpack("I", byte_string[footer_frame])[0]

        return data, footer

        # frames_in_buf = 0
        # index = 0
        # data = np.zeros(self.max_samples)  # 20 samples packed in a frame

        # while index < self.max_samples:

        #     dataFrame = self.stream.read(self.frame_size_bytes)
        #     if not dataFrame:
        #         break
        #         # if not self.open_next():
        #         #     break

        #         # dataFrame = self.stream.read(self.frame_size_bytes)
        #         # if not dataFrame:
        #         #     break

        #     dataFooter = unpack_from("I", dataFrame, self.frame_size_bytes - 4)

        #     # Check that there are no skipped frames
        #     frameCount = dataFooter[0] & self.footer_idx_samp_mask
        #     difCount = frameCount - self.last_frame
        #     if (difCount != 1):
        #         print(
        #             f"Ch [{self.ch_id}] is missing frames at {frameCount}, difference of {difCount}, check {index}\n")
        #     self.last_frame = frameCount

        #     for ptrSamp in range(0, 60, 3):
        #         tmpSampleTupple = unpack(">i", dataFrame[ptrSamp:ptrSamp + 3] + b'\x00')
        #         data[index] = tmpSampleTupple[0] * self._scale_factor
        #         index += 1

        #     frames_in_buf += 1

        #     if self.report_hw_sat:
        #         satCount = (dataFooter[0] & self.footer_sat_mask) >> 24
        #         if satCount:
        #             print ("Ch [%s] Frame %d has %d saturations" %
        #                     (self.ch_id, frameCount, satCount))

        # return data[0:index]

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


class DecimatedSegmentedReader(_TSReaderBase):
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


class DecimatedContinuousReader(_TSReaderBase):
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
