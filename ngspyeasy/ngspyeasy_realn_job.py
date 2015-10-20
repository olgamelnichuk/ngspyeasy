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


def log_exception(ex):
    get_logger(LOGGER_NAME).exception(ex)


def fix_file_permissions(projects_home, row):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id(), get_logger(LOGGER_NAME))


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
        ngspyeasy_realn_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        log_exception(ex)
        sys.exit(1)


def ngspyeasy_realn_job(tsv_conf, projects_home, sample_id, task):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = [x for x in rows2run if x.sample_id() == sample_id]

    for row in rows2run:
        try:
            run_realn(row, projects_home, task)
        finally:
            fix_file_permissions(projects_home, row)


def run_realn(row, projects_home, task):
    log_info("Running Realign job (SAMPLE_ID='%s', REALN='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.realn(), task, row.genomebuild()))

    if row.realn() == "no-realn":
        log_info("[%s] Skipping Indel Realignment for sample: '%s'" % (row.realn(), row.sample_id()))
        return

    realn_data = sample.realn_data(row, projects_home)

    if os.path.isfile(realn_data.dupl_mark_realn_bam()):
        log_info("Skipping Indel Realignment. Looks like you already ran it: %s" % realn_data.dupl_mark_realn_bam())
        return

    genomebuild = select_genomebuild(row, projects_home)

    callables = {
        "bam-realn|no-task": bam_realn,
        "gatk-realn|no-task": gatk_realn
    }

    callables.get(row.realn() + "|" + task, unrecognized_options)(row, task, realn_data, genomebuild)


def unrecognized_options(row, task, *args):
    if callable is None:
        raise ValueError(
            "Unrecognised REALN options (REALN='%s', TASK='%s')" % (row.realn(), task))


def select_genomebuild(row, projects_home):
    genomebuild = genome_build.select(row.genomebuild(), projects_home)
    if genomebuild.known_indels() is None:
        raise ValueError(
            "No genome selected for sample '%s'. GENOMEBUILD value is '%s', but it should be one of [b37, hg19]" %
            row.genomebuild())
    return genomebuild


def bam_realn(row, task, realn_data, genomebuild):
    log_debug("bam_realn (SAMPLE_ID='%s', REALN='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.realn(), task, row.genomebuild()))

    run_script("bam-realn", "bam-realn.tmpl.sh", **common_script_params(realn_data, genomebuild))


def gatk_realn(row, task, realn_data, genomebuild):
    log_debug("gatk_realn (SAMPLE_ID='%s', REALN='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.realn(), task, row.genomebuild()))

    params = dict(
        DUPEMARK_BAM_FOR_INDER_REALN_INTERVALS=realn_data.reports_path(
            realn_data.bam_prefix() + ".dupemk.bam.ForIndelRealigner.intervals")
    )

    params.update(**common_script_params(realn_data, genomebuild))

    run_script("gatk-realn", "gatk-realn.tmpl.sh", **params)


def common_script_params(realn_data, genomebuild):
    if not os.path.isfile(realn_data.dupl_mark_bam()):
        raise IOError(
            "Can't proceed with Indel Realignment as input bam doesn't exist: %s" % realn_data.dupl_mark_bam())

    row = realn_data.row()
    return dict(
        NCPU=str(row.ncpu()),
        DUPEMARK_BED=realn_data.dupl_mark_bed(),
        DUPEMARK_BAM=realn_data.dupl_mark_bam(),
        CHROMS="${chroms}",
        REFFASTA=genomebuild.ref_fasta(),
        KNOWN_INDELS=genomebuild.known_indels(),
        DUPEMARK_REALN_BAM=realn_data.dupl_mark_realn_bam(),
        DUPEMARK_REALN_FLAGSTAT=realn_data.dupl_mark_realn_bam_flagstat(),
        DUPEMARK_REALN_BED=realn_data.dupl_mark_realn_bed(),
        TMP_DIR=realn_data.tmp_dir()
    )


def run_script(dir, scriptname, **kwargs):
    base_dir = os.path.dirname(__file__)
    template_path = os.path.join(base_dir, "resources", "realn", dir, scriptname)

    log_debug("Using script template file: %s" % template_path)
    log_debug("Script params: %s" % kwargs)
    script = script_from_template(template_path)
    script.add_variables(**kwargs)
    run_command(script.to_temporary_file(), get_logger(LOGGER_NAME))


if __name__ == '__main__':
    main(sys.argv[1:])
