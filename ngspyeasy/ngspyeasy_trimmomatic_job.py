#!/usr/bin/env python

import getopt
import subprocess
import sys
import tsv_config
import projects_dir

from cmdline_options import check_cmdline_options
from logger import init_logger, log_error, log_info, log_debug


def usage():
    print """
Usage:  ngspyeasy_trimmomatic_job -c <config_file> -d <project_directory> -i <sample_id>

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
        ngspyeasy_trimmomatic_job(tsv_conf, projects_home, sample_id)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_trimmomatic_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_trimmomatic(row, projects_home)


def run_trimmomatic(row, projects_home):
    if (row.trim() == "atrim"):
        run_trimmomatic_atrim(row, projects_home)

    # else if (row.trim() == "btrim"):
    #    run_trimmomatic_btrim(row, projects_home)

    else:
        raise ValueError("Unknown trimmomatic type: %s" % row.trim())


def run_trimmomatic_atrim(row, projects_home):
    if row.genomebuild() == "b37":
        adapter_fa = "/home/pipeman/ngs_projects/ngseasy_resources/reference_genomes_b37/contaminant_list.fa"
    elif row.genomebuild() == "hg19":
        adapter_fa = "/home/pipeman/ngs_projects/ngseasy_resources/reference_genomes_hg19/contaminant_list.fa"
    elif row.genomebuild() == "hs37d5":
        adapter_fa = "/home/pipeman/ngs_projects/ngseasy_resources/reference_genomes_hs37d5/contaminant_list.fa"
    elif row.genomebuild() == "hs38DH":
        adapter_fa = "/home/pipeman/ngs_projects/ngseasy_resources/reference_genomes_hs38DH/contaminant_list.fa"
    else:
        raise ValueError("Unknown GENOMEBUILD value: '%s'" % row.genomebuild())

    fastq = [] # TODO
    log_info("Running Trimmomatic tool...")
    cmd = ["java", "-XX:ParallelGCThreads=1", "-jar", "/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar",
           "PE", "-threads", row.ncpu()] + fastq + [
              "ILLUMINACLIP:" + adapter_fa + ":2:30:10:5:true",
              "LEADING:3",
              "TRAILING:3",
              "SLIDINGWINDOW:4:15",
              "AVGQUAL:2",
              "MINLEN:75"]

    proc = subprocess.Popen(["/bin/bash", "-i", "-c", "source ~/.bashrc; " + " ".join(cmd)],
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
