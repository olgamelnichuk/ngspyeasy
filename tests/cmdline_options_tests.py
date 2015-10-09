#!/usr/bin/env python

import os
import unittest
import tempfile
from ngspyeasy.cmdargs import check_cmdline_options


class CommandLineOptionsTest(unittest.TestCase):
    def test_no_options(self):
        self.assertError(check_cmdline_options(None, None))

    def test_projects_dir_is_none(self):
        self.assertError(check_cmdline_options("config.tsv", None))

    def test_tsv_name_is_none(self):
        self.assertError(check_cmdline_options(None, "/some_dir"))

    def test_projects_dir_not_exists(self):
        self.assertError(check_cmdline_options("config.tsv", "/some_dir"))

    def test_tsv_config_not_exists(self):
        self.assertError(check_cmdline_options("config.tsv", tempfile.gettempdir()))

    def test_all_options(self):
        ngs_projects_dir = os.path.join(tempfile.gettempdir(), "ngs_projects")
        ngs_projects_config_dir = os.path.join(ngs_projects_dir, "config_files")
        if not os.path.exists(ngs_projects_config_dir):
            os.makedirs(ngs_projects_config_dir)

        file = tempfile.NamedTemporaryFile(dir=ngs_projects_config_dir)
        (tsv_name, projects_home, errmsg) = check_cmdline_options(file.name, ngs_projects_dir)
        self.assertEqual(tsv_name, os.path.basename(file.name))
        self.assertEqual(projects_home, ngs_projects_dir)
        self.assertTrue(errmsg is None)

        file = tempfile.NamedTemporaryFile(dir=ngs_projects_config_dir)
        (tsv_name, projects_home, errmsg) = check_cmdline_options(os.path.basename(file.name), ngs_projects_dir)
        self.assertEqual(tsv_name, os.path.basename(file.name))
        self.assertEqual(projects_home, ngs_projects_dir)
        self.assertTrue(errmsg is None)

    def assertError(self, (tsv_name, projects_home, errmsg)):
        self.assertTrue(tsv_name is None)
        self.assertTrue(projects_home is None)
        self.assertTrue(errmsg is not None)
