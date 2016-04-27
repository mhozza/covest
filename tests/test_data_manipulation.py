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
    def assertListEqualApprox(self, a, b, threshold):
        if len(a) > len(b):
            b += [0]*(len(a) - len(b))
        if len(a) < len(b):
            a += [0]*(len(b) - len(a))
        for i, j in zip(a, b):
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
        max_hist = max(histogram.keys())
        true_hist = [0, 2949548, 255275, 108173, 36802, 10079, 2303, 452, 79, 12, 3, 1]
        hist_l = [histogram.get(i, 0) for i in range(max_hist+1)]
        sampled_hist = sample_hist(hist_l, 2)
        self.assertListEqualApprox(sampled_hist, true_hist, 2)

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
        max_hist = max(histogram.keys())
        true_hist = [0, 2936731, 229577, 82298, 23544, 7025, 1946, 279, 23, 1]
        hist_l = [histogram.get(i, 0) for i in range(max_hist+1)]
        sampled_hist = sample_hist(hist_l, 2)
        self.assertListEqualApprox(sampled_hist, true_hist, 2)


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
        max_hist = max(histogram.keys())
        hist_l = [histogram.get(i, 0) for i in range(max_hist + 1)]
        trimmed_hist, tail =trim_hist(hist_l, 10)
        self.assertEqual(len(trimmed_hist), 10)
        self.assertEqual(sum(trimmed_hist) + tail, sum(hist_l))

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
        max_hist = max(histogram.keys())
        true_hist = [0, 2936731, 229577, 82298, 23544, 7025, 1946, 279, 23, 1]
        hist_l = [histogram.get(i, 0) for i in range(max_hist + 1)]
        trimmed_hist, tail = trim_hist(hist_l, 10)
        self.assertEqual(len(trimmed_hist), 9)
        self.assertEqual(sum(trimmed_hist) + tail, sum(hist_l))


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
