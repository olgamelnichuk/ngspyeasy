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
import job_scheduler
import tsv_config
from logger import init_logger, get_logger


def main(argv):
    parser = argparse.ArgumentParser(description="Python version of NGSeasy pipeline.")
    parser.add_argument("init", nargs="?")
    parser.add_argument("fastqc", nargs="?")
    parser.add_argument("trimmomatic", nargs="?")
    parser.add_argument("alignment", nargs="?")
    parser.add_argument("realign", nargs="?")
    parser.add_argument("bsqr", nargs="?")
    parser.add_argument("variant-calling", nargs="?")

    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory_path,
                        help="ngs_projects directory path")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")

    args = parser.parse_args(argv)

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
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
            scheduler = job_scheduler.JobScheduler()
            scheduler.start()
        except Exception, e:
            logger.exception(e)
            sys.exit(1)

    retcode = 0
    dependencies = dict()
    verbose = args.verbose
    try:
        if args.init:
            ngspyeasy_init(tsv_conf, projects_home)
        elif args.fastqc:
            ngspyeasy_fastqc(tsv_conf, projects_home, dependencies, verbose)
        elif args.trimmomatic:
            ngspyeasy_trimmomatic(tsv_conf, projects_home, dependencies, verbose)
        elif args.alignment:
            ngspyeasy_alignment(tsv_conf, projects_home, dependencies, verbose)
        elif args.realign:
            ngspyeasy_realn(tsv_conf, projects_home, dependencies, verbose)
        elif args.bsqr:
            ngspyeasy_bsqr(tsv_conf, projects_home, dependencies, verbose)
        elif args.variant_calling:
            ngspyeasy_variant_calling(tsv_conf, projects_home, dependencies, verbose)
        else:
            ngspyeasy(tsv_conf, projects_home, dependencies, verbose)
    except Exception as e:
        logger.exception(e)
        job_scheduler.stop()
        retcode = 1

    while True:
        threads = threading.enumerate()
        if len(threads) == 1:
            break
        for t in threads:
            if t != threading.currentThread():
                t.join(1)

    logger.info("Exit(retcode=%d)", retcode)
    sys.exit(retcode)


def log_info(msg):
    get_logger().info(msg)


def log_debug(msg):
    get_logger().debug(msg)


def ngspyeasy(tsv_conf, projects_home, dependencies, verbose):
    ngspyeasy_init(tsv_conf, projects_home)

    ngspyeasy_fastqc(tsv_conf, projects_home, dependencies, verbose)
    ngspyeasy_trimmomatic(tsv_conf, projects_home, dependencies, verbose)
    ngspyeasy_alignment(tsv_conf, projects_home, dependencies, verbose)
    ngspyeasy_realn(tsv_conf, projects_home, dependencies, verbose)
    ngspyeasy_bsqr(tsv_conf, projects_home, dependencies, verbose)
    ngspyeasy_variant_calling(tsv_conf, projects_home, dependencies, verbose)


def ngspyeasy_init(tsv_conf, projects_home):
    log_info("Initiating project...")
    projects_home.init_structure(tsv_conf, get_logger())

    log_info("Checking raw FastQ files...")
    projects_home.check_fastq(tsv_conf, get_logger())


def ngspyeasy_fastqc(tsv_conf, projects_home, dependencies, verbose):
    log_info("Submitting FastQC jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        cmd = JobCommand("ngspyeasy_fastqc_job.py", tsv_conf.filename(), sample_id, verbose=verbose)
        submit(cmd, "fastqc", "compbio/ngseasy-fastqc", projects_home, dependencies)


def ngspyeasy_trimmomatic(tsv_conf, projects_home, dependencies, verbose):
    log_info("Submitting Trimmomatic jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        trim_type = row.trim()

        if trim_type == "no-trim":
            log_info("[%s] No trimmomatic jobs to be run for sample: '%s'. NOT RECOMMENDED" % (trim_type, sample_id))
            continue

        if trim_type not in ["atrim", "btrim"]:
            raise ValueError("Unknown trimmomatic type: %s" % trim_type)

        cmd = JobCommand("ngspyeasy_trimmomatic_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        submit(cmd.with_task("trimmomatic"), "trim", "compbio/ngseasy-trimmomatic",
               projects_home, dependencies)

        submit(cmd.with_task("fastqc"), "trim_fastqc", "compbio/ngseasy-fastqc",
               projects_home, dependencies)


def ngspyeasy_alignment(tsv_conf, projects_home, dependencies, verbose):
    log_info("Submitting Alignment jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        aligner_type = row.aligner()

        cmd = JobCommand("ngspyeasy_alignment_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        if aligner_type == "no-align":
            log_info("[%s] No alignment jobs to be run for sample: '%s'." % (aligner_type, sample_id))
            continue
        elif aligner_type == "bwa":
            submit(cmd, "bwa", "compbio/ngseasy-bwa", projects_home, dependencies)
        elif aligner_type == "novoalign":
            submit(cmd, "novalign", "compbio/ngseasy-novoalign", projects_home, dependencies)
        elif aligner_type == "stampy":
            submit(cmd.with_task("bwa"), "stampy_bwa", "compbio/ngseasy-bwa", projects_home, dependencies)
            submit(cmd.with_task("stampy"), "stampy_stampy", "compbio/ngseasy-bwa", projects_home, dependencies)
            submit(cmd.with_task("picard_"), "stampy_picard_cleansam", "compbio/ngseasy-picardtools", projects_home,
                   dependencies)
            submit(cmd.with_task("picard_"), "stampy_picard_addorreplacereadgroups", "compbio/ngseasy-picardtools",
                   projects_home, dependencies)
        elif aligner_type == "bowtie2":
            submit(cmd, "bowtie2", "compbio/ngseasy-bowtie2", projects_home, dependencies)
        elif aligner_type == "snap":
            submit(cmd, "snap", "compbio/ngseasy-snap", projects_home, dependencies)
        else:
            raise ValueError("Unknown aligner type: %s" % aligner_type)


def ngspyeasy_realn(tsv_conf, projects_home, dependencies, verbose):
    log_info("Submitting realignment jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        realn_type = row.realn()

        cmd = JobCommand("ngspyeasy_realn_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        if realn_type == "no-realn":
            log_info("[%s] Skipping Indel Realignment for sample: '%s'." % (realn_type, sample_id))
            continue
        elif realn_type == "bam-realn":
            submit(cmd, "bam_realn", "compbio/ngseasy-glia", projects_home, dependencies)
        elif realn_type == "gatk-realn":
            submit(cmd, "gtalk_realn", "compbio/ngseasy-gatk", projects_home, dependencies)
        else:
            raise ValueError("Unknown realign type: %s" % realn_type)


def ngspyeasy_bsqr(tsv_conf, projects_home, dependencies, verbose):
    log_info("Submitting base quality score recalibration jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        bsqr_type = row.bsqr()

        cmd = JobCommand("ngspyeasy_bsqr_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        if bsqr_type == "no-bsqr":
            log_info("[%s] Skipping Base quality score recalibration for sample: '%s'" % (bsqr_type, sample_id))
            continue
        elif bsqr_type == "bam-bsqr":
            submit(cmd, "bam_bsqr", "compbio/ngseasy-base", projects_home, dependencies)
        elif bsqr_type == "gatk-bsqr":
            submit(cmd, "gatk_bsqr", "compbio/ngseasy-gatk", projects_home, dependencies)
        else:
            raise ValueError("Unknown bsqr type: %s" % bsqr_type)


def ngspyeasy_variant_calling(tsv_conf, projects_home, dependencies, verbose):
    log_info("Submitting Variant Calling jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        vc_type = row.varcaller()

        cmd = JobCommand("ngspyeasy_vc_job.py", tsv_conf.filename(), sample_id, verbose=verbose)

        submit(cmd.with_task("vc_init"), "vc_init", "compbio/ngseasy-base", projects_home, dependencies)

        if vc_type == "freebayes-parallel":
            submit(cmd, vc_type, "compbio/ngseasy-freebayes", projects_home, dependencies)
        elif vc_type == "freebayes-default":
            submit(cmd, vc_type, "compbio/ngseasy-freebayes", projects_home, dependencies)
        elif vc_type == "platypus":
            submit(cmd, vc_type, "compbio/ngseasy-platypus", projects_home, dependencies)
        elif vc_type == "platypus-default":
            submit(cmd, vc_type, "compbio/ngseasy-platypus", projects_home, dependencies)
        elif vc_type == "UnifiedGenotyper":
            submit(cmd, vc_type, "compbio/ngseasy-gatk", projects_home, dependencies)
        elif vc_type == "HaplotypeCaller":
            submit(cmd, vc_type, "compbio/ngseasy-gatk", projects_home, dependencies)
        elif vc_type == "ensemble":
            submit(cmd.with_task("freebayes"), vc_type, "compbio/ngseasy-freebayes", projects_home, dependencies)
            submit(cmd.with_task("platypus"), vc_type, "compbio/ngseasy-platypus", projects_home, dependencies)
            submit(cmd.with_task("HaplotypeCaller"), vc_type, "compbio/ngseasy-gatk", projects_home, dependencies)
        else:
            raise ValueError("Unknown variant calling type: %s" % vc_type)


def submit(command, tag, image, projects_home, dependencies):
    image += ":%s" % NGSEASYVERSION
    job_id = job_id_generator.get_next([tag, command.sample_id])
    prev_job_ids = [x for x in [dependencies.get(command.sample_id, None)] if x is not None]

    log_debug(
        "Submit job(sample_id='%s', job_id='%s', dependencies='%s')" % (command.sample_id, job_id, prev_job_ids))

    log_debug(
        "Submit job(docker_image=%s, cmd=%s)" % (image, command.cmd()))

    job_scheduler.submit(
        job_id,
        docker_cmd(job_id, image, command.cmd(), projects_home.root(), projects_home.resources_dir(), pipeman=False),
        prev_job_ids)
    dependencies[command.sample_id] = job_id


def signal_handler(signum):
    log_info("Got SIGINT(%s) signal" % str(signum))
    job_scheduler.stop()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])


class JobCommand(object):
    def __init__(self, executable, config_name, sample_id, **kwargs):
        self.executable = executable
        self.sample_id = sample_id
        self.config_name = config_name
        self.verbose = kwargs.get("verbose", False)
        self.task = kwargs.get("task", None)

    def with_task(self, task):
        return JobCommand(self.executable, self.config_name, self.sample_id, verbose=self.verbose, task=task)

    def cmd(self, ):
        cmd = ["python /ngspyeasy/bin/%s" % self.executable, "-v" if self.verbose else "", "-c", self.config_name, "-d",
               "/home/pipeman/ngs_projects", "-i", self.sample_id]
        if self.task:
            cmd += ["-t", self.task]
        return " ".join(cmd)
