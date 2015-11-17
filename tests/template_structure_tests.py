#!/usr/bin/env python

import unittest
import itertools

from ngspyeasy import pipeline_tools

class TemplateStructureTests(unittest.TestCase):
    def test_vars(self):
        for r, n in zip(pipeline_tools.refs(), itertools.count()):
            print "[%s]" % str(n), r
            #templ = pipeline_tools.find(r)



