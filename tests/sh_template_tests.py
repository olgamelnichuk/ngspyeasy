#!/usr/bin/env python

import unittest

from ngspyeasy import sh_template


class ShTemplateTest(unittest.TestCase):
    def test_realn_bam_template(self):
        tmpl = sh_template.load("realn", "bam-realn", "bam-realn.tmpl.sh")
        vars = ["NCPU", "DUPEMARK_BED", "DUPEMARK_BAM", "REFFASTA", "KNOWN_INDELS", "DUPEMARK_REALN_BAM",
                "DUPEMARK_REALN_FLAGSTAT", "DUPEMARK_REALN_BED", "TMP_DIR"]
        params = dict(zip(vars, vars))
        tmpl.validate(**params)
        print("\n".join(tmpl.source(**params)))
