#!/usr/bin/env python

import getopt
import os
import subprocess
import sys
import re
from ngspyeasy import tsv_config
from cmdline_options import check_cmdline_options
from logger import init_job_logger, log_error, log_info
from project_structure import get_config_path, get_log_dir, get_sample_dir, get_sample_fastq_path, \
    get_sample_tmp_dir, get_sample_fastq_dir


def usage():
    print """
Usage:  ngspyeasy_fastqc_job -c <config_file> -d <project_directory> -i <sample_id>

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

    log_name = tsv_name + "_fastqc_job_" + (sample_id if sample_id is not None else "all_samples")
    init_job_logger(get_log_dir(projects_home), log_name, verbose)

    tsv_conf = tsv_config.parse(get_config_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    try:
        ngspyeasy_fastqc_job(tsv_conf, projects_home, sample_id)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_fastqc_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.get_sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_fastqc(row, projects_home)


def run_fastqc(row, projects_home):
    sample_dir = get_sample_dir(projects_home, row.get_project_id(), row.get_sample_id())
    fastq1 = get_sample_fastq_path(sample_dir, row.get_fastq1())
    fastq2 = get_sample_fastq_path(sample_dir, row.get_fastq2())

    if not os.path.isfile(fastq1):
        raise OSError("File not found: %s", fastq1)

    if not os.path.isfile(fastq2):
        raise OSError("File not found: %s", fastq2)

    prefix_fastq1, type1 = get_fastqc_basename(fastq1)
    prefix_fastq2, type2 = get_fastqc_basename(fastq2)

    if prefix_fastq1 is None:
        raise ValueError("Unsupported Fastq file format: %s", fastq1)

    if prefix_fastq2 is None:
        raise ValueError("Unsupported Fastq file format: %s", fastq2)

    if type1 != type2:
        raise ValueError("Fastqc file formans are not the same")

    if type1 == "illumina":
        fqout1 = get_sample_fastq_path(sample_dir, prefix_fastq1 + "_fastqc.html")
        fqout2 = get_sample_fastq_path(sample_dir, prefix_fastq2 + "_fastqc.html")
    else:
        fqout1 = get_sample_fastq_path(sample_dir, prefix_fastq1 + "_1_fastqc.html")
        fqout2 = get_sample_fastq_path(sample_dir, prefix_fastq2 + "_2_fastqc.html")

    log_info("Check if FastQC Data already exists: [%s] and [%s]", fqout1, fqout2)
    if os.path.isfile(fqout1) and os.path.isfile(fqout2):
        log_info("FastQC Data already exists...skipping this bit")
        return

    cmd = ["/usr/local/pipeline/FastQC/fastqc", "--threads", "2", "--extract",
           "--dir", get_sample_tmp_dir(sample_dir), "--outdir", get_sample_fastq_dir(sample_dir), fastq1, fastq2]

    try:
        output = subprocess.check_output(["/bin/bash", "-c", "source ~/.bashrc && env"],
                                         stderr=subprocess.STDOUT, shell=True,
                                         env=os.environ.copy())
        retcode = 0
    except subprocess.CalledProcessError, ex:
        log_error("Command [[\n%s\n]] failed. See logs for details", " ".join(ex.cmd))
        output = ex.output
        retcode = ex.returncode

    log_info("cmd: \n" + output)
    sys.exit(retcode)


def get_fastqc_basename(fastq_file):
    illumina_patterns = [r'(.*_L.*_R[1,2]_[0-9][0-9][0-9])\.fastq\.gz',
                         r'(.*_L.*_R[1,2]_[0-9][0-9][0-9][0-9])\.fastq\.gz']

    other_patterns = [r'(.*)_[1,2]\.fastq\.gz',
                      r'(.*)_R[1,2]\.fastq\.gz',
                      r'(.*)_[1,2]\.fq\.gz',
                      r'(.*)_R[1,2]\.fq\.gz']

    for p in illumina_patterns:
        match = re.match(p, fastq_file)
        if match:
            return match.group(1), "illumina"

    for p in other_patterns:
        match = re.match(p, fastq_file)
        if match:
            return match.group(1), "other"

    return None, "unknown"


if __name__ == '__main__':
    main(sys.argv[1:])
