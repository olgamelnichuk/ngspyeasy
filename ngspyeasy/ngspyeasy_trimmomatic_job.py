#!/usr/bin/env python
import sys

import cmdargs
from shutils import run_command
import os
import tsv_config
import projects_dir
import sample
import genome_build
from logger import init_logger, get_logger

LOGGER_NAME = "trimmomatic"


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).debug(msg)


def log_error(msg):
    get_logger(LOGGER_NAME).error(msg)


def main(argv):
    args = cmdargs.parse_job_args(argv, "Trimmomatic")

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
        ngspyeasy_trimmomatic_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_trimmomatic_job(tsv_conf, projects_home, sample_id, task):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_trimmomatic(row, projects_home, task)


def run_trimmomatic(row, projects_home, task):
    log_info("Running Trimmomatic job (SAMPLE_ID='%s', TRIM='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.trim(), row.genomebuild()))

    if row.trim() not in ["atrim", "btrim", "no-trim"]:
        raise ValueError("Unrecognised TRIM option. Should be one of [atrim] [btrim] or [no-trim]: '%s'" % row.trim())

    if row.trim() == "no-trim":
        log_info("[%s] Skipping quality control of raw fastq reads. NOT RECOMMENDED" % row.trim())
        return

    if task not in ["trimmomatic", "fastqc"]:
        raise ValueError("Unknown trimmomatic task: %s" % task)

    if task == "trimmomatic":
        run_trimmomatic_task(row, projects_home)
    if task == "fastqc":
        run_fastqc_task(row, projects_home)


def run_fastqc_task(row, projects_home):
    trim_data = sample.trimmomatic_data(row, projects_home)

    fastq_files = trim_data.paired_fastq() + trim_data.unpaired_fastq()
    not_exist = filter(lambda x: not os.path.isfile(x), fastq_files)

    if len(not_exist) != 0:
        log_info("Can't proceed with (post Trimmomatic) FastQC as fastq files do not exist: %s" % not_exist)
        return

    log_info("Running (post Trimmomatic) FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "4",
           "--extract",
           "--quiet",
           "--dir", trim_data.tmp_dir(),
           "--outdir", trim_data.fastq_dir()] + fastq_files
    run_command(cmd, get_logger(LOGGER_NAME))


def run_trimmomatic_task(row, projects_home):
    trim_data = sample.trimmomatic_data(row, projects_home)

    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome is None:
        raise ValueError("Unknown GENOMEBUILD value: '%s'" % row.genomebuild())

    adapter_fa = genome.adapter_fa()

    pe = trim_data.paired_fastq()
    ue = trim_data.unpaired_fastq()
    trimmomatic_results = [pe[0], ue[0], pe[1], ue[1]]
    log_info("Checking if Trimmomatic data already exists: %s" % trimmomatic_results)

    not_exist = filter(lambda x: not os.path.isfile(x), trimmomatic_results)
    if len(not_exist) == 0:
        log_info("Trimmomatic data already exists...skipping this bit")
        return

    log_info("Running Trimmomatic tool...")
    trimmomatic_options = ["LEADING:3",
                           "TRAILING:3",
                           "SLIDINGWINDOW:4:15",
                           "AVGQUAL:2",
                           "MINLEN:75"]

    if row.trim() == "atrim":
        log_info("TRIM set to '%s' - adaptor trim. Adaptor and read quality trimming" % row.trim())
        trimmomatic_options = ["ILLUMINACLIP:" + adapter_fa + ":2:30:10:5:true"] + trimmomatic_options

    if row.trim() == "btrim":
        log_info("TRIM set to '%s' - basic trim. Just read quality trimming. No adaptor trimming" % row.trim())

    log_info("Trimmomatic options:\n %s" % "\n".join(trimmomatic_options))

    cmd = ["java", "-XX:ParallelGCThreads=1", "-jar", "/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar",
           "PE",
           "-threads", row.ncpu()] + trim_data.fastq_files() + trimmomatic_results + trimmomatic_options
    run_command(cmd, get_logger(LOGGER_NAME))


if __name__ == '__main__':
    main(sys.argv[1:])
