#!/usr/bin/env python
import sys

import cmdargs
from shutils import run_command
import os
import sample
import projects_dir
import tsv_config
from logger import init_logger, get_logger

LOGGER_NAME = "fastqc"


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).debug(msg)


def log_error(msg):
    get_logger(LOGGER_NAME).error(msg)


def log_exception(ex):
    get_logger(LOGGER_NAME).exception(ex)


def main(argv):
    args = cmdargs.parse_job_args(argv, "FastQC")

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

    retcode = 0
    try:
        ngspyeasy_fastqc_job(tsv_conf, projects_home, args.sample_id)
    except Exception as ex:
        log_exception(ex)
        retcode = 1
    finally:
        log_info("Chmod 0775 on everything under %s" % projects_home.root())
        projects_home.chmod(0775)

    sys.exit(retcode)


def ngspyeasy_fastqc_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_fastqc(row, projects_home)


def run_fastqc(row, projects_home):
    log_info("Running FastQC job (SAMPLE_ID='%s', FASTQC='%s')" % (row.sample_id(), row.fastqc()))

    if row.fastqc() == "no-fastqc":
        log_info("[%s] Skipping FastQC..." % row.fastqc())
        return

    fq_data = sample.fastqc_data(row, projects_home)

    for fastq_file in fq_data.fastq_files():
        if not os.path.isfile(fastq_file):
            raise IOError("FastQ file not found: %s", fastq_file)

    fastqc_results = fq_data.fastqc_htmls()
    log_info("Checking if FastQC results already exist: %s" % fastqc_results)

    not_exist = filter(lambda x: not os.path.isfile(x), fastqc_results)
    if len(not_exist) == 0:
        log_info("FastQC results already exist...skipping this bit")
        return

    log_info("Running FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "2",
           "--extract",
           "--dir", fq_data.tmp_dir(),
           "--outdir", fq_data.fastq_dir()] + fq_data.fastq_files()

    run_command(cmd, get_logger(LOGGER_NAME))


if __name__ == '__main__':
    main(sys.argv[1:])
