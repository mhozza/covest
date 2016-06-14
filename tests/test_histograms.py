#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_covest
----------------------------------

Tests for `covest` module.
"""

import unittest

from covest.histogram import compute_coverage_apx, sample_histogram, trim_hist


class TestComputeCovApx(unittest.TestCase):
    def test_compute_coverage_apx(self):
        histogram = {
            1: 2909,
            2: 10891,
            3: 28824,
            4: 56698,
            5: 92099,
            6: 122998,
            7: 137748,
            8: 137507,
            9: 124723,
            10: 100866,
            11: 72467,
            12: 47639,
            13: 29893,
            14: 17119,
            15: 9026,
            16: 4713,
            17: 2077,
            18: 767,
            19: 288,
            20: 139,
            21: 49,
            22: 27,
            23: 37,
            24: 16,
        }
        c, e = compute_coverage_apx(histogram, 21, 100)
        self.assertAlmostEqual(c, 10, delta=1)
        self.assertAlmostEqual(e, 0, delta=0.01)


class TestSampleHistogram(unittest.TestCase):
    def assertDictEqualApprox(self, a, b, threshold):
        for k in set(a) | set(b):
            i = a.get(k, 0)
            j = b.get(k, 0)
            if abs(i - j) > threshold:
                raise AssertionError(
                    '{} and {} differs by more then threshold {}'.format(i, j, threshold)
                )

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
        true_hist = {
            1: 2949548,
            2: 255275,
            3: 108173,
            4: 36802,
            5: 10079,
            6: 2303,
            7: 452,
            8: 79,
            9: 12,
            10: 3,
            11: 1,
        }
        sampled_hist = sample_histogram(histogram, 2)
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
        true_hist = {
            1: 2936731, 2: 229577, 3: 82298, 4: 23544, 5: 7025, 6: 1946, 7: 279, 8: 23, 9: 1
        }
        sampled_hist = sample_histogram(histogram, 2)
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
        trimmed_hist, tail = trim_hist(histogram, 10)
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
