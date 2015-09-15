#!/usr/bin/env python

import unittest
import os
from os.path import isdir
from ngspyeasy import tsv_config
from ngspyeasy.project_structure import get_config_path, get_project_dir, get_sample_dir, get_sample_tmp_dir, \
    get_sample_config_dir, get_sample_log_dir, get_sample_fastq_dir, get_sample_alignments_dir, get_sample_reports_dir, \
    get_sample_vcf_dir
from testsettings import get_resource_path
from ngspyeasy.project_structure import init_project, init_fastq
from testutils import create_ngs_projects_dir, file_permissions


class NGSProjectsStructureTest(unittest.TestCase):
    def test_directory_structure(self):
        projects_home = create_ngs_projects_dir("test_directory_structure",
                                                get_resource_path("ngspyeasy_test.config.tsv"))
        tsv_conf = tsv_config.parse(get_config_path(projects_home, "ngspyeasy_test.config.tsv"))
        init_project(tsv_conf, projects_home)

        row0 = tsv_conf.row_at(0)
        self.assertTrue(isdir(get_project_dir(projects_home, row0.get_project_id())))

        sample_dir = get_sample_dir(projects_home, row0.get_project_id(), row0.get_sample_id())
        self.assertTrue(isdir(sample_dir))
        self.assertTrue(isdir(get_sample_tmp_dir(sample_dir)))
        self.assertTrue(isdir(get_sample_config_dir(sample_dir)))
        self.assertTrue(isdir(get_sample_log_dir(sample_dir)))
        self.assertTrue(isdir(get_sample_fastq_dir(sample_dir)))
        self.assertTrue(isdir(get_sample_alignments_dir(sample_dir)))
        self.assertTrue(isdir(get_sample_reports_dir(sample_dir)))
        self.assertTrue(isdir(get_sample_vcf_dir(sample_dir)))

        for root, dirs, files in os.walk(projects_home):
            for name in files:
                self.assertTrue(0775, file_permissions(os.path.join(root, name)))

            for name in dirs:
                self.assertTrue(0775, file_permissions(os.path.join(root, name)))

    def test_fastq_not_exist(self):
        projects_home = create_ngs_projects_dir("test_fastq_not_exist",
                                                get_resource_path("ngspyeasy_test.config.tsv"))
        tsv_conf = tsv_config.parse(get_config_path(projects_home, "ngspyeasy_test.config.tsv"))
        init_project(tsv_conf, projects_home)

        try:
            init_fastq(tsv_conf, projects_home)
            self.fail("Exception should be raised, as fastq files do not exist")
        except ValueError:
            pass

    def test_fastq_exist(self):
        projects_home = create_ngs_projects_dir("test_fastq_exist",
                                                get_resource_path("ngspyeasy_test.config.tsv"))
        tsv_conf = tsv_config.parse(get_config_path(projects_home, "ngspyeasy_test.config.tsv"))
        init_project(tsv_conf, projects_home)

        raw_fastq_dir = os.path.join(projects_home, "raw_fastq")
        os.makedirs(raw_fastq_dir)
        open(os.path.join(raw_fastq_dir, "illumina.100bp.pe.wex.150x_1.fastq.gz"),'a').close()
        open(os.path.join(raw_fastq_dir, "illumina.100bp.pe.wex.150x_2.fastq.gz"),'a').close()

        try:
            init_fastq(tsv_conf, projects_home)
        except ValueError:
            self.fail("Exception should not be raised, as fastq files exist")