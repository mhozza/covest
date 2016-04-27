#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_data_manipulation
----------------------------------

Tests for `covest.data` module.
"""

import unittest
from pathlib import Path

from covest.data import load_histogram, InvalidFormatException, sample_hist, trim_hist


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


class TestSampleHistogram(unittest.TestCase):
    def assertDictEqualApprox(self, a, b, threshold):
        for k in set(a) | set(b):
            i = a.get(k, 0)
            j = b.get(k, 0)
            if abs(i - j) > threshold:
                raise AssertionError('{} and {} differs by more then threshold {}'.format(i, j, threshold))

    def test_sample_histogram(self):
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
        true_hist = {1: 2949548, 2: 255275, 3: 108173, 4: 36802, 5: 10079, 6: 2303, 7: 452, 8: 79, 9: 12, 10: 3, 11: 1}
        sampled_hist = sample_hist(histogram, 2)
        self.assertDictEqualApprox(sampled_hist, true_hist, 2)

    def test_sample_sparse_histogram(self):
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
        true_hist = {1: 2936731, 2: 229577, 3: 82298, 4: 23544, 5: 7025, 6: 1946, 7: 279, 8: 23, 9: 1}
        sampled_hist = sample_hist(histogram, 2)
        self.assertDictEqualApprox(sampled_hist, true_hist, 2)


class TestTrimHistogram(unittest.TestCase):
    def test_trim_histogram(self):
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
        trimmed_hist, tail =trim_hist(histogram, 10)
        self.assertEqual(max(trimmed_hist), 9)
        self.assertEqual(sum(trimmed_hist.values()) + tail, sum(histogram.values()))

    def test_trim_sparse_histogram(self):
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
        trimmed_hist, tail = trim_hist(histogram, 10)
        self.assertEqual(max(trimmed_hist), 8)
        self.assertEqual(sum(trimmed_hist.values()) + tail, sum(histogram.values()))


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
