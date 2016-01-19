#!/usr/bin/env python

import unittest
from ngspyeasy.utils import uniq_set


class NGSProjectsStructureTest(unittest.TestCase):
    def test_uniq_set(self):
        test1 = uniq_set(['a', 'a', 'a'])
        self.assertEqual({'a'}, set(test1))

        test2 = uniq_set(['a', 'b', 'a'])
        self.assertEqual({'a', 'b'}, set(test2))
