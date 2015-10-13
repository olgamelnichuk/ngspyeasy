#!/usr/bin/env python
import argparse
import sys

from ngspyeasy import cmdargs
from shutils import run_command
import os
import sample_data
import projects_dir
import tsv_config
from logger import init_logger, get_logger

LOGGER_NAME = "fastqc_job"


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).debug(msg)


def log_error(msg):
    get_logger(LOGGER_NAME).error(msg)


def main(argv):
    parser = argparse.ArgumentParser(description="FastQC Job")
    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory_path,
                        help="ngs_projects directory path")
    parser.add_argument("-i", "--sample_id", dest="sample_id", help="sample_id to run job on")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")

    args = parser.parse_args(argv)

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    init_logger(log_file, args.verbose, LOGGER_NAME)
    log_info("Starting...")
    log_debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    log_info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        log_error(e)
        sys.exit(1)

    try:
        ngspyeasy_fastqc_job(tsv_conf, projects_home, args.sample_id)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_fastqc_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_fastqc(row, projects_home)


def run_fastqc(row, projects_home):
    log_info("Running FastQC job (SAMPLE_ID='%s', FASTQC='%s')" % (row.sample_id(), row.fastqc()))

    sample = sample_data.create(row, projects_home)

    fastqc_results = sample.fastqc_results()
    log_info("Checking if FastQC results already exist: %s" % fastqc_results)

    not_exist = filter(lambda x: not os.path.isfile(x), fastqc_results)
    if len(not_exist) == 0:
        log_info("FastQC results already exist...skipping this bit")
        return

    log_info("Running FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "2",
           "--extract",
           "--dir", sample.tmp_dir(),
           "--outdir", sample.fastq_dir()] + sample.fastq_files()

    run_command(cmd, get_logger(LOGGER_NAME))


if __name__ == '__main__':
    main(sys.argv[1:])
