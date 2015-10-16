#!/usr/bin/env python

import argparse

import os.path


def parse_job_args(argv, name=""):
    parser = argparse.ArgumentParser(description="NGSPyEasy %s Job" % name)
    parser.add_argument("-c", "--config", dest="config", required=True, type=path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=existed_directory_path,
                        help="ngs_projects directory path")
    parser.add_argument("-i", "--sample_id", dest="sample_id", help="sample_id to run %s on" % name)
    parser.add_argument("-t", "--task", dest="task", default="no-task", help="%s task")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")
    return parser.parse_args(argv)


def existed_directory_path(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError('%s is not an existed directory path' % path)
    return path


def path_basename(path):
    return os.path.basename(path)
