#!/usr/bin/env python

import unittest
import os
from os.path import isdir
from ngspyeasy import tsv_config
from ngspyeasy import projects_dir
from testsettings import get_resource_path
from testutils import create_ngs_projects_dir, file_permissions


class NGSProjectsStructureTest(unittest.TestCase):
    def test_directory_structure(self):
        ngs_projects = create_ngs_projects_dir("test_directory_structure",
                                                get_resource_path("ngspyeasy_test.config.tsv"))

        projects_home = projects_dir.ProjectsDir(ngs_projects)

        tsv_conf = tsv_config.parse(projects_home.config_path("ngspyeasy_test.config.tsv"))
        projects_home.init_structure(tsv_conf)

        row0 = tsv_conf.row_at(0)
        self.assertTrue(isdir(projects_home.project_dir(row0.project_id())))

        sample_dir = projects_dir.SampleDir(projects_home.sample_dir(row0.project_id(), row0.sample_id()))
        self.assertTrue(isdir(sample_dir.root()))
        self.assertTrue(isdir(sample_dir.tmp_dir()))
        self.assertTrue(isdir(sample_dir.config_dir()))
        self.assertTrue(isdir(sample_dir.log_dir()))
        self.assertTrue(isdir(sample_dir.fastq_dir()))
        self.assertTrue(isdir(sample_dir.alignments_dir()))
        self.assertTrue(isdir(sample_dir.reports_dir()))
        self.assertTrue(isdir(sample_dir.vcf_dir()))

        for root, dirs, files in os.walk(projects_home.root()):
            for name in files:
                self.assertTrue(0775, file_permissions(os.path.join(root, name)))

            for name in dirs:
                self.assertTrue(0775, file_permissions(os.path.join(root, name)))

    def test_fastq_not_exist(self):
        ngs_projects = create_ngs_projects_dir("test_fastq_not_exist",
                                                get_resource_path("ngspyeasy_test.config.tsv"))
        projects_home = projects_dir.ProjectsDir(ngs_projects)

        tsv_conf = tsv_config.parse(projects_home.config_path("ngspyeasy_test.config.tsv"))
        projects_home.init_structure(tsv_conf)

        try:
            projects_home.check_fastq(tsv_conf)
            self.fail("Exception should be raised, as fastq files do not exist")
        except IOError:
            pass

    def test_fastq_exist(self):
        ngs_projects = create_ngs_projects_dir("test_fastq_exist",
                                                get_resource_path("ngspyeasy_test.config.tsv"))

        projects_home = projects_dir.ProjectsDir(ngs_projects)

        tsv_conf = tsv_config.parse(projects_home.config_path("ngspyeasy_test.config.tsv"))
        projects_home.init_structure(tsv_conf)

        os.makedirs(projects_home.raw_fastq_dir())
        open(projects_home.raw_fastq_path("illumina.100bp.pe.wex.150x_1.fastq.gz"), 'a').close()
        open(projects_home.raw_fastq_path("illumina.100bp.pe.wex.150x_2.fastq.gz"), 'a').close()

        try:
            projects_home.check_fastq(tsv_conf)
        except IOError:
            self.fail("Exception should not be raised, as fastq files exist")