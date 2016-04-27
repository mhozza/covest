#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_covest
----------------------------------

Tests for `covest` module.
"""

import unittest

from covest.covest import compute_coverage_apx


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


if __name__ == '__main__':
    import sys

    sys.exit(unittest.main())
