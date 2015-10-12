#!/usr/bin/env python
import argparse
import sys

from ngspyeasy import cmdargs
from shutils import run_command
import os
import docker
import tsv_config
import projects_dir
import sample_data
from logger import init_logger, get_logger

LOGGER_NAME = "trimmomatic_job"


def main(argv):
    parser = argparse.ArgumentParser(description="Trimmomatic Job")
    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory_path,
                        help="ngs_projects directory path")
    parser.add_argument("-i", "--sample_id", dest="sample_id", help="sample_id to run trimmomatic on")
    parser.add_argument("-t", "--task", dest="task", required=True, help="trimmomatic task: [trimmomatic | fastqc]")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")

    args = parser.parse_args(argv)

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    logger = init_logger(log_file, args.verbose, LOGGER_NAME)
    logger.info("Starting...")
    logger.debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    logger.info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        logger.error(e)
        sys.exit(1)

    try:
        ngspyeasy_trimmomatic_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        logger.error(ex)
        sys.exit(1)


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).info(msg)


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
    sample = sample_data.create(row, projects_home)

    fastq_files = sample.trimmomatic_paired_results() + sample.trimmomatic_unpaired_results()
    not_exist = filter(lambda x: not os.path.isfile(x), fastq_files)

    if len(not_exist) != 0:
        log_info("Can't proceed with (post Trimmomatic) FastQC as fastq files do not exist: %s" % not_exist)
        return

    log_info("Running (post Trimmomatic) FastQC tool...")

    cmd = ["/usr/local/pipeline/FastQC/fastqc",
           "--threads", "4",
           "--extract"
           "--quiet",
           "--dir", sample.tmp_dir(),
           "--outdir", sample.fastq_dir()] + fastq_files
    run_command(cmd, get_logger(LOGGER_NAME))


def run_trimmomatic_task(row, projects_home):
    sample = sample_data.create(row, projects_home)

    if row.genomebuild() == "b37":
        adapter_fa = docker.NGS_RESOURCES + "/reference_genomes_b37/contaminant_list.fa"
    elif row.genomebuild() == "hg19":
        adapter_fa = docker.NGS_RESOURCES + "/reference_genomes_hg19/contaminant_list.fa"
    elif row.genomebuild() == "hs37d5":
        adapter_fa = docker.NGS_RESOURCES + "/reference_genomes_hs37d5/contaminant_list.fa"
    elif row.genomebuild() == "hs38DH":
        adapter_fa = docker.NGS_RESOURCES + "/reference_genomes_hs38DH/contaminant_list.fa"
    else:
        raise ValueError("Unknown GENOMEBUILD value: '%s'" % row.genomebuild())

    pe = sample.trimmomatic_paired_results()
    ue = sample.trimmomatic_unpaired_results()
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
           "-threads", row.ncpu()] + sample.fastq_files() + trimmomatic_results + trimmomatic_options
    run_command(cmd, get_logger(LOGGER_NAME))


if __name__ == '__main__':
    main(sys.argv[1:])
