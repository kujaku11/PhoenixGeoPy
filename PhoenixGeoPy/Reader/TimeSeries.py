"""Module to read and parse native Phoenix Geophysics data formats of the MTU-5C Family

This module implements Streamed readers for segmented-decimated continuus-decimated
and native sampling rate time series formats of the MTU-5C family.
"""

__author__ = 'Jorge Torres-Solis'

# =============================================================================
# Imports
# =============================================================================

from pathlib import Path
import numpy as np
from datetime import datetime

from struct import unpack_from, unpack
import string
from PhoenixGeoPy.Reader.DataScaling import DataScaling
# =============================================================================
class Header:
    def __init__(self, **kwargs):
        self.report_hw_sat = False
        self.header_size = 128
        self.ad_plus_minus_range = 5.0  # differential voltage range that the A/D can measure (Board model dependent)
        self.channel_main_gain = None     # The value of the main gain of the board
        self.intrinsic_circuitry_gain = None  # Circuitry Gain not directly configurable by the user
        self.total_circuitry_gain = None  # Total board Gain both intrinsic gain and user-seletable gain in circuit
        self.total_selectable_gain = None  # Total of the gain that is selectable by the user (i.e. att * pre * gain)
        self.lp_frequency = None                   # Nominal cutoff freq of the configured LPF of the channel
        self.preamp_gain = 1.0
        self.attenuator_gain = 1.0
        self._recording_id = None
        self._channel_id = None
        
        self._unpack_dict = {
            "file_type": {"dtype": 'B', "index": 0},
            "file_version": {"dtype": 'B', "index": 1},
            "header_length": {"dtype": 'H', "index": 2},
            "instrument_type": {"dtype": '8s', "index": 4},
            "serial_number": {"dtype": 'cccccccc', "index": 12},
            "recording_id": {"dtype": 'I', "index": 20},
            "channel_id": {"dtype": 'B', "index": 24},
            "file_sequence": {"dtype": 'I', "index": 25},
            "frag_period": {"dtype": 'H', "index": 29},
            "ch_board_model": {"dtype": '8s', "index": 31},
            "ch_board_serial": {"dtype": '8s', "index": 39},
            "ch_firmware": {"dtype": 'I', "index": 47},
            "hardware_configuration": {"dtype": 'BBBBBBBB', "index": 51},
            "sample_rate_base":{"dtype": 'H', "index": 59},
            "sample_rate_exp":{"dtype": 'b', "index": 61},
            "bytes_per_sample":{"dtype": 'B', "index": 62},
            "frame_size":{"dtype": 'I', "index": 63},
            "decimation_node_id":{"dtype": 'H', "index": 67},
            "frame_rollover_count":{"dtype": 'H', "index": 69},
            "gps_long":{"dtype": 'f', "index": 71},
            "gps_lat":{"dtype": 'f', "index": 75},
            "gps_height":{"dtype": 'f', "index": 79},
            "gps_hacc":{"dtype": 'I', "index": 83},
            "gps_vacc":{"dtype": 'I', "index": 87},
            "timing_status":{"dtype": 'BBH', "index": 91},
            "future1":{"dtype": 'b', "index": 95},
            "future2":{"dtype": 'i', "index": 97},
            "saturated_frames":{"dtype": 'H', "index": 101},
            "missing_frames":{"dtype": 'H', "index": 103},
            "battery_voltage_mv":{"dtype": 'H', "index": 105},
            "min_signal":{"dtype": 'f', "index": 107},
            "max_signal":{"dtype": 'f', "index": 111},
            }
        
        for key, value in kwargs.items():
            setattr(self, key, value)
            
    def __str__(self):
        lines = []
        for key in sorted(self.__dict__.keys()):
            lines.append(f"{key}: {getattr(self, key)}")
            
        return "\n".join(lines)
    
    def __repr__(self):
        return self.__str__()
    
    def _has_header(self):
        if self._header is not None:
            return True
        return False
    
        
    def _unpack_value(self, key):
        if self._has_header():
            return unpack_from(
                self._unpack_dict[key]["dtype"], 
                self._header,
                self._unpack_dict[key]["index"]
                )
    
    @property
    def file_type(self):
        if self._has_header():
            return self._unpack_value("file_type")[0]
    
    @property
    def file_version(self):
        if self._has_header():
            return self._unpack_value("file_version")[0]
        
    @property
    def header_length(self):
        if self._has_header():
            return self._unpack_value("header_length")[0]
    
    @property
    def instrument_type(self):
        if self._has_header():
            return self._unpack_value("instrument_type")[0].decode("utf-8").strip(' ').strip('\x00')
    
    @property
    def instrument_serial_number(self, header):
        if self._has_header():      
            return b''.join(self._unpack_value("serial_number")).strip(b'\x00')
    
    @property
    def recording_id(self):
        if self._recording_id is None:
            if self._has_header():
                return self._unpack_value("recording_id")[0]
        else:
            return self._recording_id

    @recording_id.setter
    def recording_id(self, value):
        self._recording_id = value
        
    @property
    def start_time(self):
        if self.recording_id is not None:
            return datetime.fromtimestamp(self.recording_id).isoformat()
        return None
        
    @property
    def channel_id(self):
        if self._channel_id is None:
            if self._has_header():
                return self._unpack_value("channel_id")[0]
        else:
            return self._channel_id

    @channel_id.setter
    def channel_id(self, value):
        self._channel_id = value
    
    @property
    def file_sequence(self):
        if self._has_header():
            return self._unpack_value("file_sequence")[0]

    @property
    def frag_period(self):
        if self._has_header():
            return self._unpack_value("frag_period")[0]  
        
    @property
    def ch_board_model(self):
        if self._has_header():
            return self._unpack_value("ch_board_model")[0].decode("utf-8").strip(' ')
        
    @property
    def board_model_main(self):
        if self._has_header():
            return self.ch_board_model[0:5]

    @property
    def board_model_revision(self):
        if self._has_header():
            return self.ch_board_model[6:1]    
    
    @property
    def ch_board_serial(self):
        if self._has_header():
            value = self._unpack_value("ch_board_serial")[0].decode("utf-8").strip('\x00')
            # handle the case of backend < v0.14, which puts '--------' in ch_ser
            if all(chars in string.hexdigits for chars in value):
                return int(value, 16)
            else:
                return 0 
        
    @property
    def ch_firmware(self):
        if self._has_header():
            return self._unpack_value("ch_firmware")[0]
        
    @property
    def hardware_configuration(self):
        if self._has_header():
            return self._unpack_value("hardware_configuration") 

    @property 
    def channel_type(self):
        if self._has_header():
            if self.hardware_configuration[1] & 0x08 == 0x08:
                return "E"
            else:
                return "H"
            
    @property
    def detected_channel_type(self):
        # Channel type detected by electronics
        # this normally matches self.channel_type, but used in electronics design and testing
        if self._has_header():
            if self.hardware_configuration[1] & 0x20 == 0x20:
                return 'E'
            else:
                return 'H'

    @property
    def lp_frequency(self):
        if self._has_header():
            if self.hardware_configuration[0] & 0x80 == 0x80:            # LPF on
                if self.hardware_configuration[0] & 0x03 == 0x03:
                    return 10
                elif self.hardware_configuration[0] & 0x03 == 0x02:
                    if (self.board_model_main == "BCM03" or self.board_model_main == "BCM06"):
                        return 1000
                    else:
                        return 100
                elif self.hardware_configuration[0] & 0x03 == 0x01:
                    if (self.board_model_main == "BCM03" or self.board_model_main == "BCM06"):
                        return  10000
                    else:
                        return 1000
            else:                                      # LPF off
                if (self.board_model_main == "BCM03" or self.board_model_main == "BCM06"):
                    return 17800
                else:
                    return 10000
   
    def __popuate_peamp_gain(self, config_fp):
        if self.channel_type == "?":
            raise Exception("Channel type must be set before attemting to calculate preamp gain")
        preamp_on = bool(config_fp[0] & 0x10)
        self.preamp_gain = 1.0
        if self.channel_type == "E":
            if preamp_on is True:
                if self.board_model_main == "BCM01" or self.board_model_main == "BCM03":
                    self.preamp_gain = 4.0
                    if (self.board_model_revision == "L"):
                        #Account for BCM01-L experimental prototype
                        self.preamp_gain = 8.0
                else:
                    self.preamp_gain = 8.0
                    # Acount for experimental prototype BCM05-A
                    if self.ch_hwv[0:7] == "BCM05-A":
                        self.preamp_gain = 4.0
    
    def __populate_main_gain(self, config_fp):
        # BCM05-B and BCM06 introduced different selectable gains
        new_gains = True   # we asume any newer board will have the new gain banks
        if self.board_model_main == "BCM01" or self.board_model_main == "BCM03":
            # Original style 24 KSps boards and original 96 KSps boards
            new_gains = False
        if self.ch_hwv[0:7] == "BCM05-A":
            # Acount for experimental prototype BCM05-A, which also had original gain banks
            new_gains = False
        if config_fp[0] & 0x0C == 0x00:
            self.channel_main_gain = 1.0
        elif config_fp[0] & 0x0C == 0x04:
            self.channel_main_gain = 4.0
        elif config_fp[0] & 0x0C == 0x08:
            self.channel_main_gain = 6.0
            if not new_gains:
                self.channel_main_gain = 16.0
        elif config_fp[0] & 0x0C == 0x0C:
            self.channel_main_gain = 8.0
            if not new_gains:
                self.channel_main_gain = 32.0
   
    def __handle_sensor_range(self, config_fp):
        """This function will adjust the intrinsic circuitry gain based on the
           sensor range configuration in the configuration fingerprint
           
           For this, we consider that for the Electric channel, calibration path, or H-legacy
           sensors all go through a 1/4 gain stage, and then they get a virtial x2 gain from
           Single-ended-diff before the A/D. In the case of newer sensors (differential)
           instead of a 1/4 gain stage, there is only a 1/2 gain stage
           
           Therefore, in the E,cal and legacy sensor case the circuitry gain is 1/2, while for
           newer sensors it is 1
           """
        if self.channel_type == "?":
            raise Exception("Channel type must be set before attemting to calculate preamp gain")
        self.intrinsic_circuitry_gain = 0.5   
        if self.channel_type == "H":
            if config_fp[1] & 0x01 == 0x01:
                self.intrinsic_circuitry_gain = 1.0
   
    def __populate_attenuator_gain(self, config_fp):
        self.attenuator_gain = 1.0    # Asume attenuator off
        if self.channel_type == "?":
            raise Exception("Channel type must be set before attemting to calculate preamp gain")
        attenuator_on = bool(config_fp[4] & 0x01)
        if attenuator_on and self.channel_type == "E":
            new_attenuator = True  # By default assume that we are dealing with a newer types of boards
            if self.board_model_main == "BCM01" or self.board_model_main == "BCM03":
                # Original style 24 KSps boards and original 96 KSps boards
                new_attenuator = False
            if self.cs_hwv[0:7] == "BCM05-A":
                # Acount for experimental prototype BCM05-A, which also had original gain banks
                new_attenuator = False
   
            if new_attenuator:
                self.attenuator_gain = 523.0 / 5223.0
            else:
                self.attenuator_gain = 0.1
                

   
    def unpack_header(self, stream):
        if self.header_size > 0:
            # be sure to read from the beginning of the file
            stream.seek(0)
            header = stream.read(self.header_size)
        else:
            return
     
        print("INFO: Unpacking header")
        self.file_type = unpack_from('B', header, 0)[0]
        self.file_version = unpack_from('B', header, 1)[0]
        self.length = unpack_from('H', header, 2)[0]
        self.inst_type = unpack_from('8s', header, 4)[0].decode("utf-8").strip(' ').strip('\x00')
        self.inst_serial = b''.join(unpack_from('cccccccc', header, 12)).strip(b'\x00')
        self._rec_id = unpack_from('I', header, 20)[0]
        self._ch_id = unpack_from('B', header, 24)[0]
        self.file_sequence = unpack_from('I', header, 25)[0]
        self.frag_period = unpack_from('H', header, 29)[0]
        self.ch_hwv = unpack_from('8s', header, 31)[0].decode("utf-8").strip(' ')
        self.board_model_main = self.ch_hwv[0:5]
        self.board_model_revision = self.ch_hwv[6:1]
        self.ch_ser = unpack_from('8s', header, 39)[0].decode("utf-8").strip('\x00')
        # handle the case of backend < v0.14, which puts '--------' in ch_ser
        if all(chars in string.hexdigits for chars in self.ch_ser):
            self.ch_ser = int(self.ch_ser, 16)
        else:
            self.ch_ser = 0
        self.ch_fir = hex(unpack_from('I', header, 47)[0]) 
        self.conf_fp = unpack_from('BBBBBBBB', header, 51)
        
        # Channel type
        self.__populate_channel_type(self.conf_fp)
        # Electric channel Preamp
        self.__popuate_peamp_gain(self.conf_fp)
        # LPF
        self.__populate_lpf(self.conf_fp)
        # Main Gain Stage
        self.__populate_main_gain(self.conf_fp)
        # Sensor range
        self.__handle_sensor_range(self.conf_fp)
        # Electric channel attenuator
        self.__populate_attenuator_gain(self.conf_fp)
        # Board-wide gains
        self.total_selectable_gain = self.channel_main_gain * self.preamp_gain * self.attenuator_gain
        self.total_circuitry_gain = self.total_selectable_gain * self.intrinsic_circuitry_gain

        self.sample_rate_base = unpack_from('H', header, 59)[0]
        self.sample_rate_exp = unpack_from('b', header, 61)[0]
        self.sample_rate = self.sample_rate_base
        if self.sample_rate_exp != 0:
            self.sample_rate *= pow(10, self.sample_rate_exp)
        self.bytes_per_sample = unpack_from('B', header, 62)[0]
        self.frame_size = unpack_from('I', header, 63)[0]
        self.data_footer = self.frame_size >> 24
        self.frame_size_bytes = self.frame_size & 0x0ffffff
        self.decimation_node_id = unpack_from('H', header, 67)[0]
        self.frame_rollover_count = unpack_from('H', header, 69)[0]
        self.gps_long = unpack_from('f', header, 71)[0]
        self.gps_lat = unpack_from('f', header, 75)[0]
        self.gps_height = unpack_from('f', header, 79)[0]
        self.gps_hacc = unpack_from('I', header, 83)[0]
        self.gps_vacc = unpack_from('I', header, 87)[0]
        self.timing_status = unpack_from('BBH', header, 91)
        self.timing_flags = self.timing_status[0]
        self.timing_sat_count = self.timing_status[1]
        self.timing_stability = self.timing_status[2]
        self.future1 = unpack_from('b', header, 95)[0]
        self.future2 = unpack_from('i', header, 97)[0]
        self.saturated_frames = unpack_from('H', header, 101)[0]
        if self.saturated_frames & 0x80 == 0x80:
            self.saturated_frames &= 0x7F
            self.saturated_frames <<= 4
        self.missing_frames = unpack_from('H', header, 103)[0]
        self.battery_voltage_mV = unpack_from('H', header, 105)[0]
        self.min_signal = unpack_from('f', header, 107)[0]
        self.max_signal = unpack_from('f', header, 111)[0]


class _TSReaderBase(Header):
    def __init__(self, path, num_files=1, header_size=128, report_hw_sat=False):
        self._seq = None
        super().__init__(**{"header_size":header_size, "report_hw_sat": report_hw_sat})
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
            self.stream = open(filename, 'rb')
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

    def __init__(self, path, num_files=1, scale_to=DataScaling.AD_input_volts,
                 header_size=128, last_frame=0, channel_gain=0.5, ad_plus_minus_range = 5.0,
                 channel_type="E", report_hw_sat=False):
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
        self.input_plusminus_range = self.ad_plus_minus_range / self.total_circuitry_gain

        if self.data_scaling == DataScaling.AD_in_ADunits:
            self._scale_factor = 256
        elif self.data_scaling == DataScaling.AD_input_volts:
            self._scale_factor = self.ad_plus_minus_range / (2 ** 31)
        elif self.data_scaling == DataScaling.instrument_input_volts:
            self._scale_factor = self.input_plusminus_range / (2 ** 31)
        else:
            raise LookupError("Invalid scaling requested")

        # optimization variables
        self.footer_idx_samp_mask = int('0x0fffffff', 16)
        self.footer_sat_mask = int('0x70000000', 16)

    def read_frames(self, num_frames):
        frames_in_buf = 0
        _idx_buf = 0
        _data_buf = np.empty([num_frames * 20])  # 20 samples packed in a frame

        while (frames_in_buf < num_frames):

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
            if (difCount != 1):
                print ("Ch [%s] Missing frames at %d [%d]\n" %
                       (self.ch_id, frameCount, difCount))
            self.last_frame = frameCount

            for ptrSamp in range(0, 60, 3):
                # unpack expectes 4 bytes, but the frames are only 3?
                value = unpack(">i", dataFrame[ptrSamp:ptrSamp + 3] + b'\x00')[0]
                _data_buf[_idx_buf] = value * self._scale_factor
                _idx_buf += 1

            frames_in_buf += 1

            if self.report_hw_sat:
                satCount = (dataFooter[0] & self.footer_sat_mask) >> 24
                if satCount:
                    print ("Ch [%s] Frame %d has %d saturations" %
                           (self.ch_id, frameCount, satCount))

        return _data_buf
    
    @property
    def npts_per_frame(self):
        return int((self.frame_size_bytes - 4) / 3)
    
    def read(self):
            
        # first read in whole file
        self.stream.seek(self.header_size)
        byte_string = self.stream.read()
        
        n_frames = int(len(byte_string) / self.frame_size_bytes) 
        data_slices = [slice(ii*64,(ii +1) *60 + 4 * ii) for ii in range(n_frames)]
        footer_slices = [slice((ii)*60,(ii) *60 + 4) for ii in range(1, n_frames + 1)]
        
        data = np.zeros(n_frames * self.npts_per_frame, dtype=np.int32)
        footer = np.zeros(n_frames)
        for data_frame, footer_frame, ii in zip(data_slices, footer_slices, range(n_frames)):
            # need to make this part more efficient
            data[ii * self.npts_per_frame: (ii +1) * self.npts_per_frame] = [
                unpack(">i", byte_string[data_frame][slice(jj * 3, (jj+1) * 3)] + b"\x00")[0]
                for jj in range(self.npts_per_frame)]
            footer[ii]= unpack("I", byte_string[footer_frame])[0]
            
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
        while (bytes_to_skip > 0):
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
            self.subheader['timestamp'] = 0
            self.subheader['samplesInRecord'] = 0
        else:
            self.subheader['timestamp'] = unpack_from('I', subheaderBytes, 0)[0]
            self.subheader['samplesInRecord'] = unpack_from('I', subheaderBytes, 4)[0]
            self.subheader['satCount'] = unpack_from('H', subheaderBytes, 8)[0]
            self.subheader['missCount'] = unpack_from('H', subheaderBytes, 10)[0]
            self.subheader['minVal'] = unpack_from('f', subheaderBytes, 12)[0]
            self.subheader['maxVal'] = unpack_from('f', subheaderBytes, 16)[0]
            self.subheader['avgVal'] = unpack_from('f', subheaderBytes, 20)[0]

    def read_record_data(self):
        ret_array = np.empty([0])
        if (self.stream is not None
                and self.subheader['samplesInRecord'] is not None
                and self.subheader['samplesInRecord'] != 0):
            ret_array = np.fromfile(self.stream, dtype=np.float32, count=self.subheader['samplesInRecord'])
            if ret_array.size == 0:
                if not self.open_next():
                    return np.empty([0])
                # Array below will contain the data, or will be an np.empty array if end of series as desired
                ret_array = np.fromfile(self.stream, dtype=np.float32, count=self.subheader['samplesInRecord'])

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
                ret_array = np.append(ret_array,
                                   np.fromfile(self.stream,
                                            dtype=np.float32,
                                            count=(numSamples - ret_array.size)))
        return ret_array
