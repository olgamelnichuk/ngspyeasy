#!/usr/bin/env python

import os.path


def get_log_dir(projects_home):
    return os.path.join(projects_home, "run_log")


def get_config_dir(projects_home):
    return os.path.join(projects_home, "config_files")


def get_config_path(projects_home, tsv_name):
    return os.path.join(get_config_dir(projects_home), tsv_name)


def init_project(tsv_conf, projects_home, log):
    print "TODO"



