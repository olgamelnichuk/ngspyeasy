#!/usr/bin/env python
import os
import sys
import getopt
import signal
import threading

import job_scheduler
import tsv_config
from cmdline_options import check_cmdline_options
from project_structure import get_log_dir, get_config_path
from logger import log_error, log_info, init_main_logger, log_exception
from ngspyeasy_initiate_project import ngspyeasy_initiate_project
from ngspyeasy_initiate_fastq import ngspyeasy_initiate_fastq
from ngspyeasy_fastqc import ngspyeasy_fastqc


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

    init_main_logger(get_log_dir(projects_home), tsv_name, verbose)

    tsv_conf = tsv_config.parse(get_config_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    retcode = 0
    try:
        scheduler = job_scheduler.JobScheduler(os.path.join(get_log_dir(projects_home), "job_scheduler.log"))
        scheduler.start()

        try:
            ngspyeasy(tsv_conf, projects_home)

        except Exception, e:
            log_exception(e)
            job_scheduler.stop()
            retcode = 1

        while True:
            threads = threading.enumerate()
            if len(threads) == 1: break
            for t in threads:
                if t != threading.currentThread():
                    t.join(1)

    except Exception, e:
        log_exception(e)
        retcode = 1

    log_info("Exit(%d)", retcode)
    sys.exit(retcode)


def ngspyeasy(tsv_conf, projects_home):
    ngspyeasy_initiate_project(tsv_conf, projects_home)
    ngspyeasy_initiate_fastq(tsv_conf, projects_home)
    ngspyeasy_fastqc(tsv_conf, projects_home)
    # ngspyeasy_trimmomatic.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_alignment.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_realign.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_bsqr.run_all(tsv_config, ngs_projects_dir)
    # ngspyeasy_variant_calling.run_all(tsv_config, ngs_projects_dir)


def signal_handler(signum, frame):
    job_scheduler.stop()

if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
