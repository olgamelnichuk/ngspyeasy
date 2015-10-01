#!/usr/bin/env python

import unittest
from ngspyeasy.utils import uniq_set, recognize_fastq


class NGSProjectsStructureTest(unittest.TestCase):
    def test_uniq_set(self):
        test1 = uniq_set(['a', 'a', 'a'])
        self.assertEqual({'a'}, set(test1))

        test2 = uniq_set(['a', 'b', 'a'])
        self.assertEqual({'a', 'b'}, set(test2))

    def test_fastq_naming(self):
        ret = recognize_fastq("/test/path/illumina.100bp.pe.wex.150x_1.fastq.gz")
        self.assertEqual("other", ret.type)
        self.assertEqual("/test/path/illumina.100bp.pe.wex.150x_1_fastqc.html", ret.result)