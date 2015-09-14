#!/usr/bin/env python

import sys
import getopt
import tsv_config
from cmdline_options import check_cmdline_options
from project_structure import get_log_dir, get_config_path
from logger import init_logger, log_info, log_set_current_step
import project_structure


def usage():
    print """
Usage:  ngspyeasy_initiate_project -c <config_file> -d <project_directory>

Options:
        -c  STRING  configuration file
        -d  STRING  project directory
        -v  NULL    verbose
        -h  NULL    show this message
"""


def exit_with_error(msg):
    print >> sys.stderr, "ERROR:" + msg
    sys.exit(1)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hc:d:", ["help"])
        if len(opts) == 0:
            usage()
            sys.exit(1)

    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    tsv_config_file = None
    ngs_projects_dir = None
    verbose = False
    for opt, val in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == "-c":
            tsv_config_file = val
        elif opt == "-d":
            ngs_projects_dir = val
        elif opt == "-v":
            verbose = True
        else:
            assert False, "unhandled option"

    (tsv_name, projects_home, errmsg) = check_cmdline_options(tsv_config_file, ngs_projects_dir)
    if errmsg:
        exit_with_error(errmsg)

    init_logger(get_log_dir(projects_home), tsv_name, verbose)

    tsv_conf = tsv_config.parse(get_config_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    ngspyeasy_initiate_project(tsv_conf, projects_home)


def ngspyeasy_initiate_project(tsv_config, projects_home):
    log_set_current_step("ngspyeasy_initiate_project")
    log_info("Start: Initiate Project")
    project_structure.init_project(tsv_config, projects_home)


if __name__ == '__main__':
    main(sys.argv[1:])
