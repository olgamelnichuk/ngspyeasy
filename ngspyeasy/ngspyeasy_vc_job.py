#!/usr/bin/env python
import sys

import os

import cmdargs
import projects_dir
import tsv_config
import genome_build
import sample_data
from logger import init_logger, get_logger

LOGGER_NAME = "variant_calling"


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).debug(msg)


def log_error(msg):
    get_logger(LOGGER_NAME).error(msg)


def main(argv):
    args = cmdargs.parse_job_args(argv, "Variant Calling")

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    init_logger(log_file, args.verbose, LOGGER_NAME)
    log_info("Starting...")
    log_debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    log_info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        log_error(e)
        sys.exit(1)

    try:
        ngspyeasy_vc_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_vc_job(tsv_conf, projects_home, sample_id, task):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_vc(row, projects_home, task)


def run_vc(row, projects_home, task):
    log_info("Running variant caller job (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    options = ["freebayes-parallel", "freebayes-default", "platypus", "platypus-default", "UnifiedGenotyper",
               "HaplotypeCaller", "ensemble", "no-vc"]

    if row.varcaller() not in options:
        raise ValueError(
            "Unrecognised VARCALLER option '%s'. Should be one of %s" % (row.varcaller(), options))

    if row.varcaller() == "no-vc":
        log_info("[%s] Skipping variant calling for sample: '%s'" % (row.varcaller(), row.sample_id()))
        return

    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome.known_indels() is None:
        raise ValueError(
            "No genome selected for sample '%s'. GENOMEBUILD value is '%s', but it should be one of [b37, hg19]" %
            row.genomebuild())

    log_info("Genome build selected: '%s'" % genome.refdir())

    sample = sample_data.create(row, projects_home)

    if row.realn() == "no-realn" and row.bsqr() == "no-bsqr":
        bam_file = sample.dupl_mark_bam()
        filtered_bam = sample.dupl_mark_filtered_bam()
    else:
        bam_file = sample.dupl_mark_realn_bsqr_bam()
        filtered_bam = sample.dupl_mark_realn_bsqr_filtered_bam()

    if not os.path.isfile(bam_file):
        raise IOError("Can not find [%s] for Variant Calling." % bam_file)

    mapped_reads_bed = sample.reports_path(os.path.basename(bam_file) + ".mapped.reads.bed")
    callable_regions_bed = sample.reports_path(os.path.basename(bam_file) + ".genomecov.bed")

    # TODO


if __name__ == '__main__':
    main(sys.argv[1:])
