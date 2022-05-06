# -*- coding: utf-8 -*-

from pathlib import Path

TEST_PATH = Path(__file__).parent
DATA_PATH_NATIVE = TEST_PATH.parent.joinpath("sample_data", "native")
DATA_PATH_SEGMENTED = TEST_PATH.parent.joinpath("sample_data", "segmented")

NATIVE_0 = DATA_PATH_NATIVE.joinpath("0")
NATIVE_1 = DATA_PATH_NATIVE.joinpath("1")
NATIVE_2 = DATA_PATH_NATIVE.joinpath("2")
NATIVE_4 = DATA_PATH_NATIVE.joinpath("4")

SEGMENTED_0 = DATA_PATH_SEGMENTED.joinpath("0")
SEGMENTED_1 = DATA_PATH_SEGMENTED.joinpath("1")
SEGMENTED_2 = DATA_PATH_SEGMENTED.joinpath("2")
SEGMENTED_4 = DATA_PATH_SEGMENTED.joinpath("4")

NATIVE_SEQUENCE_0 = list(NATIVE_0.glob("*.bin"))
NATIVE_SEQUENCE_1 = list(NATIVE_1.glob("*.bin"))
NATIVE_SEQUENCE_2 = list(NATIVE_2.glob("*.bin"))
NATIVE_SEQUENCE_4 = list(NATIVE_4.glob("*.bin"))
