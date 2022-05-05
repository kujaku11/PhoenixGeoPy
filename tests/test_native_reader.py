# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:11:44 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest
import numpy as np

from __init__ import NATIVE_0, NATIVE_1, NATIVE_2, NATIVE_4
from PhoenixGeoPy.readers import NativeReader

# =============================================================================


class TestNativeReader_0(unittest.TestCase):
    def setUp(self):

        self.path = NATIVE_0
        self.sequence = list(NATIVE_0.glob("*.bin"))

    def test_header(self):

        n = NativeReader(self.sequence[0])
        true_values = {
            "file_type": 1,
            "file_version": 3,
            "header_length": 128,
            "instrument_type": "MTU-5C",
            "instrument_serial_number": "10128",
            "recording_id": 1619492349,
            "channel_id": "0",
            "file_sequence": 0,
            "frag_period": 60,
            "ch_board_model": "BCM01-I",
            "ch_board_serial": 200803,
            "ch_firmware": 65567,
            "hardware_configuration": (4, 3, 0, 0, 0, 10, 128, 0),
            "sample_rate_base": 24000,
            "sample_rate_exp": 0,
            "bytes_per_sample": 3,
            "frame_size": 64,
            "decimation_node_id": 0,
            "frame_rollover_count": 0,
            "gps_long": -79.39388275146484,
            "gps_lat": 43.69647979736328,
            "gps_elevation": 97.2187728881836,
            "gps_horizontal_accuracy": 18.332,
            "gps_vertical_accuracy": 44.669,
            "timing_status": 55,
            "future1": 27,
            "future2": 0,
            "saturated_frames": 0,
            "missing_frames": 0,
            "battery_voltage_mv": 12.475,
            "min_signal": -2.054893970489502,
            "max_signal": 2.0543980598449707,
        }

        for key in true_values.keys():
            with self.subTest(key):
                true_value = true_values[key]
                read_value = getattr(n, key)
                if isinstance(true_value, (int, float)):
                    self.assertAlmostEqual(true_value, read_value)
                else:
                    self.assertEqual(true_value, read_value)


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
