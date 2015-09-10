#!/usr/bin/env python

import os.path


def projects_log_dir(projects_home):
    return os.path.join(projects_home, "run_log")


def projects_conf_dir(projects_home):
    return os.path.join(projects_home, "config_files")


def projects_conf_relpath(projects_home, tsv_name):
    return os.path.join(projects_conf_dir(projects_home), tsv_name)


def init_project(tsv_conf, projects_home, log):
    print "TODO"



