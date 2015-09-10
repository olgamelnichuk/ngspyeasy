#!/usr/bin/env python

import unittest
import tempfile
from ngspyeasy import tsv_config


class TsvConfigTest(unittest.TestCase):

    def test_empty_tsv(self):
        tsv = tempfile.TemporaryFile()
        #config = tsv_config.parse(tsv, log)
        #self.assertTrue(config.isEmpty())



