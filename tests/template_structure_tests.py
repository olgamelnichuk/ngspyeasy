#!/usr/bin/env python

import unittest
import itertools
import os

from ngspyeasy import pipeline_tools
from ngspyeasy import sh_template

class TemplateStructureTests(unittest.TestCase):
    def test_vars(self):
        for r, n in zip(pipeline_tools.refs(), itertools.count()):
            print "[%s]" % str(n), r
            tmpls = pipeline_tools.find(r)
            for t in tmpls:
                self.validate(t)

    def validate(self, t):
        tmpl = sh_template.load(os.path.join("tools", t["ref"], t["template"]))
        d = dict()
        for v in (t["vars"] + t["input"] + t["output"]):
            d[v] = "DUMMY"
        tmpl.validate(**d)




