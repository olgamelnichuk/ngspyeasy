#!/usr/bin/env python
import argparse
import sys
import signal
import threading

import cmdargs
import job_id_generator
from docker import docker_cmd
from settings import NGSEASYVERSION
import projects_dir
import os
import job_scheduler
import tsv_config
from logger import init_logger, get_logger


def main(argv):
    parser = argparse.ArgumentParser(description='Python version of NGSeasy pipeline.')
    parser.add_argument('init', nargs='?')
    parser.add_argument('fastqc', nargs='?')
    parser.add_argument('trimmomatic', nargs='?')
    parser.add_argument('alignment', nargs='?')
    parser.add_argument('realign', nargs='?')
    parser.add_argument('bsqr', nargs='?')
    parser.add_argument('variant-calling', nargs='?')

    parser.add_argument('-c', '--config', dest='config', required=True, type=cmdargs.path_basename,
                        help='TSV configuration file name')
    parser.add_argument('-d', '--projects-dir', dest='projects_dir', required=True, type=cmdargs.existed_directory_path,
                        help='ngs_projects directory path')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='turn ON verbose mode')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1', help='print software version')

    args = parser.parse_args(argv)

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.main_log_file(args.config)
    print 'Opening log file: %s' % log_file

    logger = init_logger(log_file, args.verbose)
    logger.info("Starting...")
    logger.debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    logger.info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        logger.error(e)
        sys.exit(1)

    logger.info("TSV config stats: %s" % tsv_conf.stats())
    logger.info("Starting job scheduler...")

    try:
        scheduler = job_scheduler.JobScheduler(os.path.join(projects_home.log_dir(), "job_scheduler.log"))
        scheduler.start()
    except Exception, e:
        logger.exception(e)
        sys.exit(1)

    retcode = 0
    try:
        if args.init:
            ngspyeasy_init(tsv_conf, projects_home, logger)
        elif args.fastqc:
            ngspyeasy_fastqc(tsv_conf, projects_home, logger)
        elif args.trimmomatic:
            ngspyeasy_trimmomatic(tsv_conf, projects_home, logger)
        elif args.alignment:
            ngspyeasy_alignment(tsv_conf, projects_home, logger)
        elif args.realign:
            ngspyeasy_realign(tsv_conf, projects_home, logger)
        elif args.bsqr:
            ngspyeasy_bsqr(tsv_conf, projects_home, logger)
        elif args.variant_calling:
            ngspyeasy_variant_calling(tsv_conf, projects_home, logger)
        else:
            ngspyeasy(tsv_conf, projects_home, logger)
        logger.info("All jobs have been submitted.")
    except Exception as e:
        logger.exception(e)
        job_scheduler.stop()
        retcode = 1

    while True:
        threads = threading.enumerate()
        if len(threads) == 1: break
        for t in threads:
            if t != threading.currentThread():
                t.join(1)

    logger.info("Exit(%d)", retcode)
    sys.exit(retcode)


def ngspyeasy(tsv_conf, projects_home, logger):
    ngspyeasy_init(tsv_conf, projects_home, logger)

    dependencies = dict()
    ngspyeasy_fastqc(tsv_conf, projects_home, logger, dependencies)
    ngspyeasy_trimmomatic(tsv_conf, projects_home, logger, dependencies)
    ngspyeasy_alignment(tsv_conf, projects_home, logger, dependencies)
    ngspyeasy_realign(tsv_conf, projects_home, logger, dependencies)
    ngspyeasy_bsqr(tsv_conf, projects_home, logger, dependencies)
    ngspyeasy_variant_calling(tsv_conf, projects_home, logger, dependencies)


def ngspyeasy_init(tsv_conf, projects_home, logger):
    logger.info("Initiating project...")
    projects_home.init_structure(tsv_conf, logger)

    logger.info("Checking raw FastQ files...")
    projects_home.check_fastq(tsv_conf, logger)


def ngspyeasy_fastqc(tsv_conf, projects_home, logger, dependencies=None):
    logger.info("Submitting FastQC jobs...")

    if dependencies is None:
        dependencies = dict()

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        cmd = ["python /ngspyeasy/bin/ngspyeasy_fastqc_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        submit("fastqc", sample_id, projects_home, "compbio/ngseasy-fastqc:" + NGSEASYVERSION, cmd, logger,
               dependencies)


def ngspyeasy_trimmomatic(tsv_conf, projects_home, logger, dependencies=None):
    logger.info("Submitting Trimmomatic jobs...")

    if dependencies is None:
        dependencies = dict()

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        trim_type = row.trim()

        if trim_type == "no-trim":
            logger.info("[%s] No trimmomatic jobs to be run for sample: '%s'. NOT RECOMMENDED" % (trim_type, sample_id))
            continue

        if trim_type not in ["atrim", "btrim"]:
            raise ValueError("Unknown trimmomatic type: %s" % trim_type)

        cmd = ["python /ngspyeasy/bin/ngspyeasy_trimmomatic_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        submit("trimomatic", sample_id, projects_home, "compbio/ngseasy-trimmomatic:" + NGSEASYVERSION, cmd, logger,
               dependencies)

        cmd = ["python /ngspyeasy/bin/ngspyeasy_trimmomatic_fastqc_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        submit("trimomatic_fastqc", sample_id, projects_home, "compbio/ngseasy-fastqc:" + NGSEASYVERSION, cmd, logger,
               dependencies)


def ngspyeasy_alignment(tsv_conf, projects_home, logger, dependencies=None):
    logger.ifno("Submitting Alignment jobs...")

    if dependencies is None:
        dependencies = dict()

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        aligner_type = row.aligner()

        cmd = ["python /ngspyeasy/bin/ngspyeasy_alignment_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        if aligner_type == "no-align":
            logger.info("[%s] No alignment jobs to be run for sample: '%s'." % (aligner_type, sample_id))
            continue
        elif aligner_type == "bwa":
            commands = [("compbio/ngseasy-bwa:" + NGSEASYVERSION, cmd)]
        elif aligner_type == "novoalign":
            commands = [("compbio/ngseasy-novoalign:" + NGSEASYVERSION, cmd)]
        elif aligner_type == "stampy":
            commands = [
                ("compbio/ngseasy-bwa:" + NGSEASYVERSION, cmd + ["-t", "bwa"]),
                ("compbio/ngseasy-stampy:" + NGSEASYVERSION, cmd + ["-t", "stampy"]),
                ("compbio/ngseasy-picardtools:" + NGSEASYVERSION, cmd + ["-t", "picard_cleansam"]),
                ("compbio/ngseasy-picardtools:" + NGSEASYVERSION, cmd + ["-t", "picard_addorreplacereadgroups"])
            ]
        elif aligner_type == "bowtie2":
            commands = [("compbio/ngseasy-bowtie2:" + NGSEASYVERSION, cmd)]
        elif aligner_type == "snap":
            commands = [("compbio/ngseasy-snap:" + NGSEASYVERSION, cmd)]
        else:
            raise ValueError("Unknown aligner type: %s" % aligner_type)

        for (image, command) in commands:
            submit("alignement", sample_id, projects_home, image, command, logger, dependencies)


def ngspyeasy_realign(tsv_conf, projects_home, logger, dependencies=None):
    logger.info("Submitting Realign jobs...")
    logger.info("Coming soon...")


def ngspyeasy_bsqr(tsv_conf, projects_home, logger, dependencies=None):
    logger.info("Submitting BSQR jobs...")
    logger.info("Coming soon...")


def ngspyeasy_variant_calling(tsv_conf, projects_home, logger, dependencies=None):
    logger.info("Submitting Variant Calling jobs...")
    logger.info("Coming soon...")


def submit(source, sample_id, projects_home, docker_image, cmd, logger, dependencies):
    job_id = job_id_generator.get_next([source, sample_id])
    prev_job_ids = [x for x in [dependencies.get(sample_id, None)] if x is not None]

    logger.debug(
        "Submit job(sample_id='%s', job_id='%s', dependencies='%s')" % (sample_id, job_id, prev_job_ids))

    logger.debug(
        "Submit job(docker_image=%s, cmd=%s)" % (docker_image, cmd))

    job_scheduler.submit(
        job_id, docker_cmd(job_id, docker_image, " ".join(cmd), projects_home.root(),
                           projects_home.resources_dir(), pipeman=False), prev_job_ids)
    dependencies[sample_id] = job_id


def signal_handler(signum, frame):
    get_logger().info("Got SIGINT(%s) signal" % str(signum))
    job_scheduler.stop()


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
