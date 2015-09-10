#!/usr/bin/env python

import sys
import os.path
from project_structure import projects_conf_relpath, projects_conf_dir


def check_cmdline_options(tsv_config_file, ngs_projects_dir):
    if not ngs_projects_dir:
        return None, None, "Projects directory is not specified."

    if not tsv_config_file:
        return None, None, "TSV config file is not specified."

    projects_home = os.path.abspath(ngs_projects_dir)

    if not os.path.isdir(projects_home):
        return None, None, "Projects directory '" + projects_home + "' does not exist."

    tsv_absolute_path = projects_conf_relpath(projects_home, os.path.basename(tsv_config_file))

    if os.path.abspath(tsv_config_file) != tsv_absolute_path:
        return None, None, "Config file must be in the projects config directory: '" + projects_conf_dir(
            projects_home) + "'"

    if not os.path.isfile(tsv_absolute_path):
        return None, None, "Config file '" + tsv_absolute_path + "' does not exist."

    return os.path.basename(tsv_absolute_path), projects_home, None


def exit_with_error(msg):
    print >> sys.stderr, "[ngspyeasy]:ERROR:" + msg
    sys.exit(1)
