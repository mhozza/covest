#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_data_manipulation
----------------------------------

Tests for `covest.data` module.
"""

import unittest
from pathlib import Path

from covest.data import load_histogram, InvalidFormatException


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

        hist = load_histogram(str(self.folder / 'simulated_c10_e0.05_r100_k21.hist'))
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

        hist = load_histogram(str(self.folder / 'simulated_c10_e0.05_r100_k21_sparse.hist'))
        self.assertDictEqual(hist, true_histogram)

    def test_load_nonexisting_file(self):
        with self.assertRaises(FileNotFoundError):
            load_histogram('nonexisting_histogram.hist')

    def test_load_invalid_histogram(self):
        with self.assertRaises(InvalidFormatException):
            hist = load_histogram(str(self.folder / 'invalid.hist'))

    def test_load_empty_histogram(self):
        hist = load_histogram(str(self.folder / 'empty.hist'))
        self.assertEqual(len(hist), 0)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
