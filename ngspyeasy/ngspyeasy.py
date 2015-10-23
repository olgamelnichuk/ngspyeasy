#!/usr/bin/env python
import argparse
import sys
import signal
import threading

import cmdargs
import job_id_generator
import docker
from settings import NGSEASYVERSION
import projects_dir
import job_scheduler
import tsv_config
from logger import init_logger, logger


class JobSubmitter(object):
    def __init__(self, projects_home, mode="docker", pipeman=False):
        self.dependencies = dict()
        self.projects_home = projects_home
        self.mode = mode
        self.pipeman = pipeman

    def dependencies_for(self, sample_id):
        return [x for x in [self.dependencies.get(sample_id, None)] if x is not None]

    def update_dependencies(self, sample_id, job_id):
        self.dependencies[sample_id] = job_id

    def submit(self, cmd, image, tag):
        image += ":%s" % NGSEASYVERSION
        job_id = job_id_generator.get_next([tag, cmd.sample_id])
        job_dependencies = self.dependencies_for(cmd.sample_id)

        log_debug(
            "Submit job(sample_id='%s', job_id='%s', dependencies='%s', cmd=[%s])" % (
                cmd.sample_id, job_id, job_dependencies, cmd.as_string()))

        job_scheduler.submit(
            job_id, self.wrap(job_id, job_dependencies, image, cmd.as_string()), job_dependencies)
        self.update_dependencies(cmd.sample_id, job_id)

    def wrap(self, job_id, job_dependencies, image, cmd_string):
        projects_root = self.projects_home.root()
        resources_root = self.projects_home.resources_dir()

        if self.in_lsf_mode():
            return docker.wrap_lsf(job_id, image, cmd_string, projects_root, resources_root, job_dependencies,
                                   pipeman=self.pipeman)
        return docker.wrap(job_id, image, cmd_string, projects_root, resources_root, pipeman=self.pipeman)

    def in_lsf_mode(self):
        return self.mode == "docker-lsf"


def main(argv):
    parser = argparse.ArgumentParser(description="Python version of NGSeasy pipeline.")
    parser.add_argument("init", nargs="?")
    parser.add_argument("fastqc", nargs="?")
    parser.add_argument("trimmomatic", nargs="?")
    parser.add_argument("alignment", nargs="?")
    parser.add_argument("realign", nargs="?")
    parser.add_argument("bsqr", nargs="?")
    parser.add_argument("variant_calling", nargs="?")

    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory,
                        help="ngs_projects directory path")
    parser.add_argument("-r", "--resources-dir", dest="resources_dir", type=cmdargs.existed_directory,
                        help="ngs_resources directory path")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--test", dest="test", action="store_true", help="turn ON test mode")
    # parser.add_argument("-s", "--settings", dest="settings", type=cmdargs.existed_file, help="settings file")
    parser.add_argument("-m", "--mode", dest="mode", choices=["docker", "docker-lsf"], default="docker",
                        help="running mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")

    args = parser.parse_args(argv)

    projects_home = projects_dir.ProjectsDir(args.projects_dir, args.resources_dir)
    log_file = projects_home.main_log_file(args.config)
    print "Opening log file: %s" % log_file

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

    if not args.init:
        try:
            logger.info("Starting job scheduler...")
            scheduler = job_scheduler.JobScheduler(args.test)
            scheduler.start()
        except Exception, e:
            logger.exception(e)
            sys.exit(1)

    verbose = args.verbose
    try:
        submitter = JobSubmitter(projects_home)
        if args.init:
            ngspyeasy_init(tsv_conf, projects_home)
        elif args.fastqc:
            ngspyeasy_fastqc(tsv_conf, submitter, verbose)
        elif args.trimmomatic:
            ngspyeasy_trimmomatic(tsv_conf, submitter, verbose)
        elif args.alignment:
            ngspyeasy_alignment(tsv_conf, submitter, verbose)
        elif args.realign:
            ngspyeasy_realn(tsv_conf, submitter, verbose)
        elif args.bsqr:
            ngspyeasy_bsqr(tsv_conf, submitter, verbose)
        elif args.variant_calling:
            ngspyeasy_variant_calling(tsv_conf, submitter, verbose)
        else:
            ngspyeasy(tsv_conf, projects_home, submitter, verbose)

        job_scheduler.all_done()
    except Exception as e:
        logger.exception(e)
        job_scheduler.stop()

    while True:
        threads = threading.enumerate()
        if len(threads) == 1:
            break
        for t in threads:
            if t != threading.currentThread():
                t.join(1)


def ngspyeasy(tsv_conf, projects_home, submitter, verbose):
    ngspyeasy_init(tsv_conf, projects_home)

    ngspyeasy_fastqc(tsv_conf, submitter, verbose)
    ngspyeasy_trimmomatic(tsv_conf, submitter, verbose)
    ngspyeasy_alignment(tsv_conf, submitter, verbose)
    ngspyeasy_realn(tsv_conf, submitter, verbose)
    ngspyeasy_bsqr(tsv_conf, submitter, verbose)
    ngspyeasy_variant_calling(tsv_conf, submitter, verbose)


def ngspyeasy_init(tsv_conf, projects_home):
    logger().info("Initiating project...")
    projects_home.init_structure(tsv_conf)

    logger().info("Checking raw FastQ files...")
    projects_home.check_fastq(tsv_conf)


def ngspyeasy_fastqc(tsv_conf, submitter, verbose):
    logger().info("Submitting FastQC jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        fastqc_type = row.fastqc()

        if fastqc_type not in ["fastqc", "no-fastqc"]:
            raise ValueError("Unknown fastqc type: %s" % fastqc_type)

        if fastqc_type == "no-fastqc":
            logger().info("[%s] No fastqc jobs to be run for sample: '%s'" % (fastqc_type, sample_id))
            continue

        cmd = docker.JobCommand("ngspyeasy_fastqc_job.py", tsv_conf.filename(), sample_id, verbose=verbose)
        submitter.submit(cmd, "compbio/ngseasy-fastqc", fastqc_type)


def ngspyeasy_trimmomatic(tsv_conf, submitter, verbose):
    logger().info("Submitting Trimmomatic jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        trim_type = row.trim()

        if trim_type not in ["atrim", "btrim", "no-trim"]:
            raise ValueError("Unknown trimmomatic type: %s" % trim_type)

        if trim_type == "no-trim":
            logger().info(
                "[%s] No trimmomatic jobs to be run for sample: '%s'. NOT RECOMMENDED" % (trim_type, sample_id))
            continue

        cmd = docker.JobCommand("ngspyeasy_trimmomatic_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        submitter.submit(cmd.add_task("trimmomatic"), "compbio/ngseasy-trimmomatic", trim_type)

        submitter.submit(cmd.add_task("fastqc"), "compbio/ngseasy-fastqc", trim_type)


def ngspyeasy_alignment(tsv_conf, submitter, verbose):
    logger().info("Submitting alignment jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        aligner_type = row.aligner()

        cmd = docker.JobCommand("ngspyeasy_alignment_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        if aligner_type == "no-align":
            logger().info("[%s] No alignment jobs to be run for sample: '%s'." % (aligner_type, sample_id))
            continue
        elif aligner_type == "bwa":
            submitter.submit(cmd, "compbio/ngseasy-bwa", aligner_type)
        elif aligner_type == "novoalign":
            submitter.submit(cmd, "compbio/ngseasy-novoalign", aligner_type)
        elif aligner_type == "stampy":
            submitter.submit(cmd.add_task("bwa"), "compbio/ngseasy-bwa", aligner_type)
            submitter.submit(cmd.add_task("stampy"), "compbio/ngseasy-bwa", aligner_type)
            submitter.submit(cmd.add_task("picard_cleansam"), "compbio/ngseasy-picardtools", aligner_type)
            submitter.submit(cmd.add_task("picard_addorreplacereadgroups"), "compbio/ngseasy-picardtools", aligner_type)
        elif aligner_type == "bowtie2":
            submitter.submit(cmd, "compbio/ngseasy-bowtie2", aligner_type)
        elif aligner_type == "snap":
            submitter.submit(cmd, "compbio/ngseasy-snap", aligner_type)
        else:
            raise ValueError("Unknown aligner type: %s" % aligner_type)


def ngspyeasy_realn(tsv_conf, submitter, verbose):
    logger().info("Submitting realignment jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        realn_type = row.realn()

        cmd = docker.JobCommand("ngspyeasy_realn_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        if realn_type == "no-realn":
            logger().info("[%s] Skipping Indel Realignment for sample: '%s'." % (realn_type, sample_id))
            continue
        elif realn_type == "bam-realn":
            submitter.submit(cmd, "compbio/ngseasy-glia", realn_type)
        elif realn_type == "gatk-realn":
            submitter.submit(cmd, "compbio/ngseasy-gatk", realn_type)
        else:
            raise ValueError("Unknown realign type: %s" % realn_type)


def ngspyeasy_bsqr(tsv_conf, submitter, verbose):
    logger().info("Submitting base quality score recalibration jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        bsqr_type = row.bsqr()

        cmd = docker.JobCommand("ngspyeasy_bsqr_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        if bsqr_type == "no-bsqr":
            logger().info("[%s] Skipping Base quality score recalibration for sample: '%s'" % (bsqr_type, sample_id))
            continue
        elif bsqr_type == "bam-bsqr":
            submitter.submit(cmd, "compbio/ngseasy-base", bsqr_type)
        elif bsqr_type == "gatk-bsqr":
            submitter.submit(cmd, "compbio/ngseasy-gatk", bsqr_type)
        else:
            raise ValueError("Unknown bsqr type: %s" % bsqr_type)


def ngspyeasy_variant_calling(tsv_conf, submitter, verbose, tag="varcaller"):
    logger().info("Submitting variant calling jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        vc_type = row.varcaller()

        cmd = docker.JobCommand("ngspyeasy_vc_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        submitter.submit(cmd.add_task("prepare"), "compbio/ngseasy-base", vc_type)

        if vc_type == "freebayes-parallel":
            submitter.submit(cmd, "compbio/ngseasy-freebayes", vc_type)
        elif vc_type == "freebayes-default":
            submitter.submit(cmd, "compbio/ngseasy-freebayes", vc_type)
        elif vc_type == "platypus":
            submitter.submit(cmd, "compbio/ngseasy-platypus", vc_type)
        elif vc_type == "platypus-default":
            submitter.submit(cmd, "compbio/ngseasy-platypus", vc_type)
        elif vc_type == "UnifiedGenotyper":
            submitter.submit(cmd, "compbio/ngseasy-gatk", vc_type)
        elif vc_type == "HaplotypeCaller":
            submitter.submit(cmd, "compbio/ngseasy-gatk", vc_type)
        elif vc_type == "ensemble":
            submitter.submit(cmd.add_task("freebayes"), "compbio/ngseasy-freebayes", vc_type)
            submitter.submit(cmd.add_task("platypus"), "compbio/ngseasy-platypus", vc_type)
            submitter.submit(cmd.add_task("HaplotypeCaller"), "compbio/ngseasy-gatk", vc_type)
        else:
            raise ValueError("Unknown variant calling type: %s" % vc_type)


def signal_handler(signum, frame):
    logger().info("Got SIGINT(%s) signal" % str(signum))
    job_scheduler.stop()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
