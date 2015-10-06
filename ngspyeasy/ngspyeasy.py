#!/usr/bin/env python
import sys
import getopt
import signal
import threading
import projects_dir

import os
import job_scheduler
import tsv_config
from cmdline_options import check_cmdline_options
from logger import log_info, log_exception, init_logger
from ngspyeasy_initiate_project import ngspyeasy_initiate_project
from ngspyeasy_initiate_fastq import ngspyeasy_initiate_fastq
from ngspyeasy_fastqc import ngspyeasy_fastqc
from ngspyeasy_trimmomatic import ngspyeasy_trimmomatic
from ngspyeasy_alignment import ngspyeasy_alignment


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
    print >> sys.stderr, "ERROR: " + msg
    print >> sys.stderr, "Exiting..."
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
        scheduler = job_scheduler.JobScheduler(os.path.join(projects_dir.log_dir(projects_home), "job_scheduler.log"))
        scheduler.start()
    except Exception, e:
        log_exception(e)
        sys.exit(1)

    retcode = 0
    try:
        ngspyeasy(tsv_conf, projects_home)
    except Exception, e:
        log_exception(e)
        job_scheduler.stop()
        retcode = 1

    log_info("All jobs have been submitted.")

    while True:
        threads = threading.enumerate()
        if len(threads) == 1: break
        for t in threads:
            if t != threading.currentThread():
                t.join(1)

    log_info("Exit(%d)", retcode)
    sys.exit(retcode)


def ngspyeasy(tsv_conf, projects_home):
    ngspyeasy_initiate_project(tsv_conf, projects_home)
    ngspyeasy_initiate_fastq(tsv_conf, projects_home)

    dependencies = dict()
    ngspyeasy_fastqc(tsv_conf, projects_home, dependencies)
    ngspyeasy_trimmomatic(tsv_conf, projects_home, dependencies)
    ngspyeasy_alignment(tsv_conf, projects_home, dependencies)
    # ngspyeasy_realign.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_bsqr.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_variant_calling.run_all(tsv_config, ngs_projects_dir)


def signal_handler(signum, frame):
    log_info("Got SIGINT(%s) signal" % str(signum))
    job_scheduler.stop()


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
