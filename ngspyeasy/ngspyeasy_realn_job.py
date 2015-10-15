#!/usr/bin/env python
import sys

import genome_build
import sample
from shutils import script_from_template, run_command
import os
import cmdargs
import projects_dir
import tsv_config
from logger import init_logger, get_logger

LOGGER_NAME = "realn"


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).debug(msg)


def log_error(msg):
    get_logger(LOGGER_NAME).error(msg)


def main(argv):
    args = cmdargs.parse_job_args(argv, "Indel Realignment")

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
        ngspyeasy_realn_job(tsv_conf, projects_home, args.sample_id)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_realn_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_realn(row, projects_home)


def run_realn(row, projects_home):
    log_info("Running Realign job (SAMPLE_ID='%s', REALN='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.realn(), row.genomebuild()))

    if row.realn() not in ["gatk-realn", "bam-realn", "no-realn"]:
        raise ValueError(
            "Unrecognised REALN option. Should be one of [bam-realn] [gatk-realn] or [no-realn]: '%s'" % row.realn())

    if row.realn() == "no-realn":
        log_info("[%s] Skipping Indel Realignment for sample: '%s'" % (row.realn(), row.sample_id()))
        return

    realn_data = sample.realn_data(row, projects_home)

    if os.path.isfile(realn_data.dupl_mark_realn_bam()):
        log_info("Skipping Indel Realignment. Looks like you already ran it: %s" % realn_data.dupl_mark_realn_bam())
        return

    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome.known_indels() is None:
        raise ValueError(
            "No genome selected for sample '%s'. GENOMEBUILD value is '%s', but it should be one of [b37, hg19]" %
            row.genomebuild())

    log_info("Genome build selected: '%s'" % genome.refdir())

    base_dir = os.path.dirname(__file__)
    template_path = os.path.join(base_dir, "resources", "realn", row.realn(), row.realn() + ".tmpl.sh")

    log_debug("Using script template file: %s" % template_path)

    script = script_from_template(template_path)

    log_debug("Script template to run: %s" % script.source())

    bam_prefix = realn_data.bam_prefix()
    dupl_mark_bed = realn_data.dupl_mark_bed()
    dupl_mark_bam = realn_data.dupl_mark_bam()

    if not os.path.isfile(dupl_mark_bed):
        raise IOError("Can't find mark duplication bed file: %s " % dupl_mark_bed)

    if not os.path.isfile(dupl_mark_bed):
        raise IOError("Can't find mark duplication bam file: %s " % dupl_mark_bam)

    script.add_variables(
        NCPU=str(row.ncpu()),
        DUPEMARK_BED=dupl_mark_bed,
        DUPEMARK_BAM=dupl_mark_bam,
        CHROMS="${chroms}",
        REFFASTA=genome.ref_fasta(),
        KNOWN_INDELS=genome.known_indels(),
        DUPEMARK_REALN_BAM=realn_data.dupl_mark_realn_bam(),
        DUPEMARK_REALN_FLAGSTAT=realn_data.dupl_mark_realn_bam_flagstat(),
        DUPEMARK_REALN_BED=realn_data.dupl_mark_realn_bed(),
        TMP_DIR=realn_data.tmp_dir()
    )

    if row.realn() == "gatk-realn":
        script.add_variables(
            DUPEMARK_BAM_FOR_INDER_REALN_INTERVALS=realn_data.reports_path(
                bam_prefix + ".dupemk.bam.ForIndelRealigner.intervals")
        )

    log_debug("Script template variables:\n %s" % "\n".join(script.variable_assignments()))

    run_command(script.to_temporary_file(), get_logger(LOGGER_NAME))


if __name__ == '__main__':
    main(sys.argv[1:])
