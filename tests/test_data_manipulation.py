#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_data_manipulation
----------------------------------

Tests for `covest.data` module.
"""

import os
import tempfile
import unittest

from covest.data import InvalidFormatException, load_histogram, save_histogram
from pathlib import Path


class TestLoadHistogram(unittest.TestCase):
    def setUp(self):
        self.folder = Path(__file__).parent / 'data'

    def test_load_histogram(self):
        true_histogram = {
            1: 5311035,
            2: 310620,
            3: 223796,
            4: 150045,
            5: 81833,
            6: 37443,
            7: 14628,
            8: 5075,
            9: 1505,
            10: 386,
            11: 133,
            12: 22,
            13: 5,
            14: 9,
            15: 1,
        }

        hist, meta = load_histogram(str(self.folder / 'simulated_c10_e0.05_r100_k21.hist'))
        self.assertDictEqual(hist, true_histogram)

    def test_load_sparse_histogram(self):
        true_histogram = {
            1: 5311035,
            2: 310620,
            3: 223796,
            4: 150045,
            6: 37443,
            7: 14628,
            8: 5075,
            12: 22,
            13: 5,
            15: 1,
        }

        hist, meta = load_histogram(str(self.folder / 'simulated_c10_e0.05_r100_k21_sparse.hist'))
        self.assertDictEqual(hist, true_histogram)

    def test_load_nonexisting_file(self):
        with self.assertRaises(FileNotFoundError):
            load_histogram('nonexisting_histogram.hist')

    def test_load_invalid_histogram(self):
        with self.assertRaises(InvalidFormatException):
            load_histogram(str(self.folder / 'invalid.hist'))

    def test_load_empty_histogram(self):
        hist, meta = load_histogram(str(self.folder / 'empty.hist'))
        self.assertEqual(len(hist), 0)


class TestSaveHistogram(unittest.TestCase):
    def setUp(self):
        handle, self.fname = tempfile.mkstemp()
        os.close(handle)

    def setTearDown(self):
        os.unlink(self.fname)

    def test_save_histogram(self):
        histogram = {
            1: 5311035,
            2: 310620,
            3: 223796,
            4: 150045,
            5: 81833,
            6: 37443,
            7: 14628,
            8: 5075,
            9: 1505,
            10: 386,
            11: 133,
            12: 22,
            13: 5,
            14: 9,
            15: 1,
        }

        orig_meta = {
            'test_key': 'test_value',
        }

        save_histogram(histogram, self.fname, meta=orig_meta)
        hist, meta = load_histogram(self.fname)
        self.assertDictEqual(hist, histogram)
        self.assertDictEqual(meta, orig_meta)

    def test_save_sparse_histogram(self):
        histogram = {
            1: 5311035,
            2: 310620,
            3: 223796,
            4: 150045,
            6: 37443,
            7: 14628,
            8: 5075,
            12: 22,
            13: 5,
            15: 1,
        }

        orig_meta = {
            'test_key': 'test_value',
        }

        save_histogram(histogram, self.fname, meta=orig_meta)
        hist, meta = load_histogram(self.fname)
        self.assertDictEqual(hist, histogram)
        self.assertDictEqual(meta, orig_meta)

    def test_save_no_meta(self):
        histogram = {
            1: 5311035,
            2: 310620,
            3: 223796,
            4: 150045,
            6: 37443,
            7: 14628,
            8: 5075,
            12: 22,
            13: 5,
            15: 1,
        }

        save_histogram(histogram, self.fname)
        hist, meta = load_histogram(self.fname)
        self.assertDictEqual(hist, histogram)
        self.assertDictEqual(meta, {})

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
