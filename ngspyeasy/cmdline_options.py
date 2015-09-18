#!/usr/bin/env python

import sys
import os.path

from project_structure import get_config_path, get_config_dir


def check_cmdline_options(tsv_config_file, ngs_projects_dir):
    (projects_home, errmsg) = check_ngs_projects_dir_option(ngs_projects_dir)
    if errmsg is not None:
        return None, None, errmsg

    (tsv_name, errmsg) = check_tsv_config_file_option(tsv_config_file, projects_home)
    if errmsg is not None:
        return None, None, errmsg

    return tsv_name, projects_home, None


def check_ngs_projects_dir_option(ngs_projects_dir):
    if not ngs_projects_dir:
        return None, "Projects directory is a required parameter. Can't find one."

    projects_home = os.path.abspath(ngs_projects_dir)

    if not os.path.isdir(projects_home):
        return None, "Projects directory '" + projects_home + "' does not exist."

    return projects_home, None


def check_tsv_config_file_option(tsv_config_file, projects_home):
    if not tsv_config_file:
        return None, None, "TSV config file is not specified."

    expected_path = get_config_path(projects_home, os.path.basename(tsv_config_file))

    if os.path.isabs(tsv_config_file):
        if os.path.abspath(tsv_config_file) != expected_path:
            return None, "Config file must be in the projects config directory: '" + get_config_dir(
                projects_home) + "'"

    if not os.path.isfile(expected_path):
        return None, "Config file '" + expected_path + "' does not exist."

    return os.path.basename(expected_path), None
