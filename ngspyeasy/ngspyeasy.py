#!/usr/bin/env python

import sys
import getopt
import tsv_config
from cmdline_options import check_cmdline_options
from project_structure import get_log_dir, get_config_path
from logger import init_logger, log_error

from ngspyeasy_initiate_project import ngspyeasy_initiate_project
from ngspyeasy_initiate_fastq import ngspyeasy_initiate_fastq
#from ngspyeasy_fastqc import ngspyeasy_fastqc


def usage():
    print """
Usage:  ngspyeasy -c <config_file> -d <project_directory>

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

    try:
        ngspyeasy(tsv_conf, projects_home)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy(tsv_conf, projects_home):
    ngspyeasy_initiate_project(tsv_conf, projects_home)
    ngspyeasy_initiate_fastq(tsv_conf, projects_home)
    # ngspyeasy_fastqc(tsv_conf, projects_home)
    # ngspyeasy_trimmomatic.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_alignment.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_realign.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_bsqr.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_variant_calling.run_all(tsv_config, ngs_projects_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
