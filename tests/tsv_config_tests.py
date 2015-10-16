#!/usr/bin/env python
import os.path
import unittest
import tempfile
from testsettings import get_resource_path
from ngspyeasy import tsv_config


class TsvConfigTest(unittest.TestCase):
    def test_empty_tsv(self):
        tmp = tempfile.NamedTemporaryFile()
        config = tsv_config.parse(tmp.name)

        self.assertTrue(config.is_empty())
        self.assertEqual(os.path.basename(tmp.name), config.filename())

    def test_columns(self):
        config = tsv_config.parse(get_resource_path("ngspyeasy_test.config.tsv"))

        self.assertFalse(config.is_empty())
        self.assertEqual(1, config.row_size())
        self.assertEqual(23, config.col_size())
        self.assertTrue(config.has_header())
        row = config.row_at(0)
        self.assertEqual("GCAT_Data", row.project_id())
        self.assertEqual("NA12878", row.sample_id())
        self.assertEqual("illumina.100bp.pe.wex.150x_1.fastq.gz", row.fastq1())
        self.assertEqual("illumina.100bp.pe.wex.150x_2.fastq.gz", row.fastq2())
        self.assertEqual("no-fastqc", row.fastqc())
        self.assertEqual("no-trim", row.trim())
        self.assertEqual("bwa", row.aligner())
        self.assertEqual("platypus", row.varcaller())
        self.assertEqual("WEX", row.ngs_type())
        self.assertEqual("ILLUMINA", row.ngs_platform())
        self.assertEqual("100bp150x.PE", row.dna_prep_library_id())
        self.assertEqual("hg19", row.genomebuild())
        self.assertEqual("ngspyeasy_test.config.tsv", config.filename())
