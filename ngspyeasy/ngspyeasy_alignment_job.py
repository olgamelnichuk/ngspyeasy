#!/usr/bin/env python
import getopt
import subprocess
import sys
import sample_data

import projects_dir
import tsv_config
from cmdline_options import check_cmdline_options
from logger import init_logger, log_error, log_set_current_step, log_info


def usage():
    print """
Usage:  ngspyeasy_alignment_job -c <config_file> -d <project_directory> -i <sample_id> -t <task>

Options:
        -c  STRING  configuration file
        -d  STRING  project directory
        -v  NULL    verbose
        -h  NULL    show this message
        -i  STRING sample id
        -t  STRING task name (e.g. for stampy aligner tasks are: ['stampy_bwa', 'stampy_stampy', 'stampy_picard1', 'stampy_picard2'])
"""


def exit_with_error(msg):
    print >> sys.stderr, "ERROR:" + msg
    sys.exit(1)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hvc:d:i:t:", ["help"])
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
    task = None
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
        elif opt == "-t":
            task = val
        else:
            assert False, "unhandled option"

    (tsv_name, projects_home, errmsg) = check_cmdline_options(tsv_config_file, ngs_projects_dir)
    if errmsg:
        exit_with_error(errmsg)

    init_logger(projects_dir.sample_log_file(projects_home, tsv_name, sample_id), verbose)
    log_set_current_step("ngspyeasy_fastqc_job")

    tsv_conf = tsv_config.parse(projects_dir.config_full_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    try:
        ngspyeasy_alignment_job(tsv_conf, projects_home, sample_id, task)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_alignment_job(tsv_conf, projects_home, sample_id, task=None):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_alignment(row, projects_home, task)


def run_alignment(row, projects_home, task):
    log_info("Running alignmnent job (SAMPLE_ID='%s', ALIGNER='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), row.genomebuild()))

    sample = sample_data.create(row, projects_home)

    log_info("FastQ files: %s" % sample.fastq_files())

    platform_unit = find_platform_unit(sample.fastq_files()[0])

    log_info("platform_unit='%s'" % platform_unit)


def find_platform_unit(fastq_file):
    p1 = subprocess.Popen(["zcat %s" % fastq_file], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["head -1"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(["sed 's/:/\\t/' - "], stdin=p2.stdout, stdout=subprocess.PIPE)
    p4 = subprocess.Popen(["cut -f 1"], stdin=p3.stdout, stdout=subprocess.PIPE)
    p5 = subprocess.Popen(["sed 's/@//g' - "], stdin=p4.stdout, stdout=subprocess.PIPE)

    return p5.communicate()[0]


if __name__ == '__main__':
    main(sys.argv[1:])
