#!/usr/bin/env python

import sys
import getopt
import tsv_config
from cmdline_options import check_cmdline_options
from logger import init_logger, log_info, log_error, log_set_current_step
import projects_dir


def usage():
    print """
Usage:  ngspyeasy_initiate_fastq -c <config_file> -d <project_directory>

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

    init_logger(projects_dir.main_log_file(projects_home, tsv_name), verbose)

    tsv_conf = tsv_config.parse(projects_dir.config_full_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    try:
        ngspyeasy_initiate_fastq(tsv_conf, projects_home)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_initiate_fastq(tsv_conf, projects_home):
    log_set_current_step("ngspyeasy_initiate_fastq")
    log_info("Start: Initiate Fastq")

    projects_dir.init_fastq(tsv_conf, projects_home)


if __name__ == '__main__':
    main(sys.argv[1:])