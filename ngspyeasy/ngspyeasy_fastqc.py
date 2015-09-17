#!/usr/bin/env python
import sys
import getopt
from ngspyeasy import job_scheduler
from ngspyeasy.docker import docker_cmd
from ngspyeasy.settings import NGSEASYVERSION
import tsv_config
from cmdline_options import check_cmdline_options
from project_structure import get_log_dir, get_config_path, get_resources_dir
from logger import init_logger, log_info, log_error, log_set_current_step
import job_id_generator


def usage():
    print """
Usage:  ngspyeasy_fastqc -c <config_file> -d <project_directory>

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
        ngspyeasy_fastqc(tsv_conf, projects_home)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_fastqc(tsv_conf, projects_home):
    log_set_current_step("ngspyeasy_fastqc")
    log_info("Start: FastQC")

    for row in tsv_conf.all_rows():
        sample_id = row.get_sample_id()
        cmd = ["/ngspyeasy/bin/ngspyeasy_fastqc_job", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        job_id = job_id_generator.get_next(["fastqc", sample_id])

        job_scheduler.submit(
            job_id, docker_cmd(job_id, "compbio/ngseasy-fastqc:" + NGSEASYVERSION, " ".join(cmd), projects_home,
                               get_resources_dir(projects_home), pipeman=False), [])

    if __name__ == '__main__':
        main(sys.argv[1:])
