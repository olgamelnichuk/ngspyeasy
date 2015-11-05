#!/usr/bin/env python
import sys

import cmdargs
from shcmd import run_command
import os
import sample
import projects_dir
import tsv_config
from logger import init_logger, logger


def fix_file_permissions(projects_home, row, uid, gid):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id(), uid, gid)


def main(argv):
    args = cmdargs.parse_job_args(argv, "FastQC")

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
        ngspyeasy_fastqc_job(tsv_conf, projects_home, args.sample_id, args.uid, args.gid)
    except Exception as ex:
        logger().exception(ex)
        sys.exit(1)


def ngspyeasy_fastqc_job(tsv_conf, projects_home, sample_id, uid, gid):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        try:
            run_fastqc(row, projects_home)
        finally:
            fix_file_permissions(projects_home, row, uid, gid)


def run_fastqc(row, projects_home):
    logger().info("Running FastQC job (SAMPLE_ID='%s', FASTQC='%s')" % (row.sample_id(), row.fastqc()))

    if row.fastqc() == "no-fastqc":
        logger().info("[%s] Skipping FastQC..." % row.fastqc())
        return

    if row.fastqc() != "qc-fastqc":
        raise ValueError("Unrecognised FastQC option: %s" % row.fastqc())

    fq_data = sample.fastqc_data(row, projects_home)

    for fastq_file in fq_data.fastq_files():
        if not os.path.isfile(fastq_file):
            raise IOError("FastQ file not found: %s", fastq_file)

    fastqc_results = fq_data.fastqc_htmls()
    logger().info("Checking if FastQC results already exist: %s" % fastqc_results)

    not_exist = filter(lambda x: not os.path.isfile(x), fastqc_results)
    if len(not_exist) == 0:
        logger().info("FastQC results already exist...skipping this bit")
        return

    logger().info("Running FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "2",
           "--extract",
           "--dir", fq_data.tmp_dir(),
           "--outdir", fq_data.fastq_dir()] + fq_data.fastq_files()

    run_command(cmd)


if __name__ == '__main__':
    main(sys.argv[1:])
