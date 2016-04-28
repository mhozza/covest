#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest

from covest.models import select_model, BasicModel, RepeatsModel


class TestCovest(unittest.TestCase):
    def test_model_selection_fulltext(self):
        self.assertEqual(select_model('basic'), BasicModel)
        self.assertEqual(select_model('repeat'), RepeatsModel)

    def test_model_selection_part(self):
        self.assertEqual(select_model('b'), BasicModel)
        self.assertEqual(select_model('r'), RepeatsModel)

    def test_model_selection_invalid(self):
        with self.assertRaises(ValueError):
            select_model('someInvalidModel')

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
