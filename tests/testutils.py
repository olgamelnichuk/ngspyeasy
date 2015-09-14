#!/usr/bin/env python

import os
import tempfile
from shutil import copyfile, rmtree


def create_ngs_projects_dir(name, tsv_config_file):
    ngs_projects_dir = os.path.join(tempfile.gettempdir(), name)
    if os.path.isdir(ngs_projects_dir):
        rmtree(ngs_projects_dir)

    ngs_projects_config_dir = os.path.join(ngs_projects_dir, "config_files")

    os.makedirs(ngs_projects_config_dir)

    copyfile(tsv_config_file, os.path.join(ngs_projects_config_dir, os.path.basename(tsv_config_file)))
    return ngs_projects_dir


def file_permissions(file):
    return oct(os.stat(file).st_mode & 0777)
