#!/usr/bin/env python

from mock import patch
import unittest
import tempfile
from testsettings import get_resource_path
from ngspyeasy import tsv_config


class TsvConfigTest(unittest.TestCase):

    @patch("ngspyeasy.tsv_config.log_warn")
    def test_empty_tsv(self, mock_log_warn):
        tmp = tempfile.NamedTemporaryFile()
        config = tsv_config.parse(tmp.name)

        self.assertTrue(config.is_empty())
        self.assertTrue(mock_log_warn.called)

    @patch("ngspyeasy.tsv_config.log_error")
    def test_columns(self, mock_log_error):
        config = tsv_config.parse(get_resource_path("ngspyeasy_test.config.tsv"))

        self.assertFalse(config.is_empty())
        self.assertEqual(1, config.row_size())
        self.assertEqual(23, config.col_size())
        self.assertTrue(config.has_header())
        self.assertEqual("GCAT_Data", config.row_at(0).get_project_id())
        self.assertEqual("NA12878", config.row_at(0).get_sample_id())
        self.assertEqual("illumina.100bp.pe.wex.150x_1.fastq.gz", config.row_at(0).get_fastq1())
        self.assertEqual("illumina.100bp.pe.wex.150x_2.fastq.gz", config.row_at(0).get_fastq2())
        self.assertFalse(mock_log_error.called)
