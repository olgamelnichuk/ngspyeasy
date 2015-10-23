#!/usr/bin/env python
import sys

import cmdargs
from shcmd import run_command
import os
import tsv_config
import projects_dir
import sample
import genome_build
from logger import init_logger, logger


def fix_file_permissions(projects_home, row):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id())


def main(argv):
    args = cmdargs.parse_job_args(argv, "Trimmomatic")

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    init_logger(log_file, args.verbose)
    logger().info("Starting...")
    logger().debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    logger().info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        logger().error(e)
        sys.exit(1)

    try:
        ngspyeasy_trimmomatic_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        logger().exception(ex)
        sys.exit(1)


def ngspyeasy_trimmomatic_job(tsv_conf, projects_home, sample_id, task):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        try:
            run_trimmomatic(row, projects_home, task)
        finally:
            fix_file_permissions(projects_home, row)


def run_trimmomatic(row, projects_home, task):
    logger().info("Running Trimmomatic job (SAMPLE_ID='%s', TRIM='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.trim(), task, row.genomebuild()))

    if row.trim() == "no-trim":
        logger().info("[%s] Skipping quality control of raw fastq reads. NOT RECOMMENDED" % row.trim())
        return

    callables = {
        "atrim|trimmomatic": atrim,
        "atrim|fastqc": fastqc,
        "btrim|trimmomatic": btrim,
        "btrim|fastqc": fastqc,
    }

    callables.get(row.trim() + "|" + task, unrecognized_options)(row, projects_home, task)


def unrecognized_options(row, projects_home, task):
    if callable is None:
        raise ValueError(
            "Unrecognised options (TRIM='%s', TASK='%s')" % (row.trim(), task))


def fastqc(row, projects_home, task):
    logger().debug("fastqc (SAMPLE_ID='%s', TRIM='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.trim(), task, row.genomebuild()))

    trim_data = sample.trimmomatic_data(row, projects_home)
    fastqc_reports = trim_data.trim_fastqc_htmls()
    not_exist = [x for x in fastqc_reports if not os.path.isfile(x)]
    if len(not_exist) == 0:
        logger().info("FastQC results already exist. Skipping this part... %s" % fastqc_reports)
        return

    fastq_files = trim_data.paired_fastq() + trim_data.unpaired_fastq()
    not_exist = [x for x in fastq_files if not os.path.isfile(x)]

    if len(not_exist) != 0:
        raise IOError("Can't proceed with (post Trimmomatic) FastQC as fastq files do not exist: %s" % not_exist)

    logger().info("Running (post Trimmomatic) FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "4",
           "--extract",
           "--quiet",
           "--dir", trim_data.tmp_dir(),
           "--outdir", trim_data.fastq_dir()] + fastq_files
    run_command(cmd)


def atrim(row, projects_home, task):
    logger().debug("atrim (SAMPLE_ID='%s', TRIM='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.trim(), task, row.genomebuild()))

    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome is None:
        raise ValueError("Unknown GENOMEBUILD value: '%s'" % row.genomebuild())

    adapter_fa = genome.adapter_fa()

    trim(row, projects_home, [
        "LEADING:3",
        "TRAILING:3",
        "SLIDINGWINDOW:4:15",
        "AVGQUAL:2",
        "MINLEN:75",
        "ILLUMINACLIP:" + adapter_fa + ":2:30:10:5:true"

    ])


def btrim(row, projects_home, task):
    logger().debug("btrim (SAMPLE_ID='%s', TRIM='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.trim(), task, row.genomebuild()))

    trim(row, projects_home, [
        "LEADING:3",
        "TRAILING:3",
        "SLIDINGWINDOW:4:15",
        "AVGQUAL:2",
        "MINLEN:75"
    ])


def trim(row, projects_home, options):
    trim_data = sample.trimmomatic_data(row, projects_home)

    pe = trim_data.paired_fastq()
    ue = trim_data.unpaired_fastq()
    trim_results = [pe[0], ue[0], pe[1], ue[1]]

    not_exist = [x for x in trim_results if not os.path.isfile(x)]
    if len(not_exist) == 0:
        logger().info("Trimmomatic data already exists...skipping this bit: %s" % trim_results)
        return

    logger().info("Running Trimmomatic tool...")
    logger().info("Trimmomatic options:\n %s" % "\n".join(options))

    cmd = ["java", "-XX:ParallelGCThreads=1", "-jar", "/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar",
           "PE",
           "-threads", row.ncpu()] + trim_data.fastq_files() + trim_results + options
    run_command(cmd)


if __name__ == '__main__':
    main(sys.argv[1:])
