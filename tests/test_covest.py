#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_covest
----------------------------------

Tests for `covest` module.
"""

import unittest

from covest.covest import select_model
from covest.models import BasicModel, RepeatsModel


class TestCovest(unittest.TestCase):
    def test_model_selection_fulltext(self):
        self.assertEqual(select_model('simple'), BasicModel)
        self.assertEqual(select_model('repeat'), RepeatsModel)

    def test_model_selection_part(self):
        self.assertEqual(select_model('s'), BasicModel)
        self.assertEqual(select_model('r'), RepeatsModel)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
