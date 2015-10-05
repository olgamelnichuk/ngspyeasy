#!/usr/bin/env python

import getopt
import os
import sys

import sample_data
import projects_dir
import tsv_config
from cmdline_options import check_cmdline_options, run_command
from logger import init_logger, log_error, log_info, log_set_current_step


def usage():
    print """
Usage:  ngspyeasy_trimmomatic_fastqc_job -c <config_file> -d <project_directory> -i <sample_id>

Options:
        -c  STRING  configuration file
        -d  STRING  project directory
        -v  NULL    verbose
        -h  NULL    show this message
        -i  STRING sample id
"""


def exit_with_error(msg):
    print >> sys.stderr, "ERROR:" + msg
    sys.exit(1)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hvc:d:i:", ["help"])
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
    sample_id = None
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
        elif opt == "-i":
            sample_id = val
        else:
            assert False, "unhandled option"

    (tsv_name, projects_home, errmsg) = check_cmdline_options(tsv_config_file, ngs_projects_dir)
    if errmsg:
        exit_with_error(errmsg)

    init_logger(projects_dir.sample_log_file(projects_home, tsv_name, sample_id), verbose)
    log_set_current_step("ngspyeasy_trimmomatic_fastqc_job")

    tsv_conf = tsv_config.parse(projects_dir.config_full_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    try:
        ngspyeasy_trimmomatic_fastqc_job(tsv_conf, projects_home, sample_id)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_trimmomatic_fastqc_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_fastqc(row, projects_home)


def run_fastqc(row, projects_home):
    log_info("(post Trimmomatic) FastQC Job (SAMPLE_ID='%s')" % row.sample_id())

    if not row.trim() in ["atrim", "btrim", "no-trim"]:
        raise ValueError("Unrecognised TRIM option. Should be one of [atrim] [btrim] or [no-trim]: '%s'" % row.trim())

    if row.trim() == "no-trim":
        log_info("Skipping (post Trimmomatic) FastQC...")
        return

    sample = sample_data.create(row, projects_home)
    fastq_files = sample.trimmomatic_paired_results() + sample.trimmomatic_unpaired_results()
    not_exist = filter(lambda x: not os.path.isfile(x), fastq_files)
    if len(not_exist) != 0:
        log_info("Can't proceed with (post Trimmomatic) FastQC as fastq files do not exist: %s", not_exist)
        return

    log_info("Running (post Trimmomatic) FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "4",
           "--extract"
           "--quiet",
           "--dir", sample.tmp_dir(),
           "--outdir", sample.fastq_dir()] + \
           fastq_files

    run_command(cmd)


if __name__ == '__main__':
    main(sys.argv[1:])
