# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:27:51 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import pytest
import numpy as np

from . import (
    NATIVE_SEQUENCE_0,
    NATIVE_SEQUENCE_1,
    NATIVE_SEQUENCE_2,
    NATIVE_SEQUENCE_4,
)
from PhoenixGeoPy.readers import TSReaderBase

# =============================================================================


@pytest.fixture(
    params=NATIVE_SEQUENCE_0
    + NATIVE_SEQUENCE_1
    + NATIVE_SEQUENCE_2
    + NATIVE_SEQUENCE_4,
)
def fn(request):
    yield request.param


@pytest.fixture
def base_reader(fn):
    return TSReaderBase(fn)


def test_base_path(fn, base_reader):
    assert base_reader.base_path == fn


def test_base_dir(fn, base_reader):
    assert base_reader.base_dir == fn.parent


def test_file_name(fn, base_reader):
    assert base_reader.file_name == fn.name


def test_file_extension(fn, base_reader):
    assert base_reader.file_extension == fn.suffix


def test_instrument_id(fn, base_reader):
    assert base_reader.instrument_id == fn.stem.split("_")[0]


def test_open_file(base_reader):
    assert base_reader._open_file(base_reader.base_path) == True
