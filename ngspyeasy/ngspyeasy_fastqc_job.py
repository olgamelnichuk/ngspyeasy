#!/usr/bin/env python

import getopt
import os
import subprocess
import sys

import utils
import projects_dir
import tsv_config
from cmdline_options import check_cmdline_options
from logger import init_logger, log_error, log_info, log_debug


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

    init_logger(projects_dir.sample_log_file(projects_home, tsv_name, sample_id), verbose)

    tsv_conf = tsv_config.parse(projects_dir.config_full_path(projects_home, tsv_name))
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
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_fastqc(row, projects_home)


def run_fastqc(row, projects_home):
    sample_dir = projects_dir.sample_dir(projects_home, row.project_id(), row.sample_id())
    fastq = [row.fastq1(), row.fastq2()]
    fastq = map(lambda x: projects_dir.fastq_full_path(sample_dir, x), fastq)

    for fq_file in fastq:
        if not os.path.isfile(fq_file):
            raise IOError("File does not exist: %s", fq_file)

    fastq_parsed = map(lambda x: utils.recognize_fastq(x), fastq)
    fastq_types = set(map(lambda x: x.type, fastq_parsed))

    if len(fastq_types) > 1:
        raise ValueError("Fastqc file formats are not the same: %s" % str(fastq_types))

    fastq_results = map(lambda x: x.result, fastq_parsed)
    log_info("Checking if FastQC results already exists: %s", fastq_results)

    fastq_results = filter(lambda x: not os.path.isfile(x), fastq_results)
    if len(fastq_results) == 0:
        log_info("FastQC results already exists...skipping this bit")
        return

    log_info("Running FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc", "--threads", "2", "--extract",
           "--dir", projects_dir.sample_tmp_dir(sample_dir), "--outdir",
           projects_dir.sample_fastq_dir(sample_dir)] + fastq

    proc = subprocess.Popen(["/bin/bash", "-i", "-c", "source ~/.bashrc; echo $CLASSPATH; " + " ".join(cmd)],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)

    stdout = []
    for line in iter(proc.stdout.readline, ''):
        sys.stdout.write(line)
        sys.stdout.flush()
        stdout.append(line)

    if proc.returncode:
        log_error("Command [[\n%s\n]] failed. See logs for details", " ".join(cmd))

    log_debug("cmd: \n" + "".join(stdout))
    sys.exit(proc.returncode)


if __name__ == '__main__':
    main(sys.argv[1:])
