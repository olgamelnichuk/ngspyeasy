#!/usr/bin/env python
import sys

from shcmd import run_command
import sh_template
import os
import cmdargs
import projects_dir
import tsv_config
import sample
import genome_build
from logger import init_logger, log_info, log_debug, log_error, log_exception


def fix_file_permissions(projects_home, row):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id())


def main(argv):
    args = cmdargs.parse_job_args(argv, "Base quality score recalibration")

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    init_logger(log_file, args.verbose)
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
        ngspyeasy_bsqr_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        log_exception(ex)
        sys.exit(1)


def ngspyeasy_bsqr_job(tsv_conf, projects_home, sample_id, task):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = [x for x in rows2run if x.sample_id() == sample_id]

    for row in rows2run:
        try:
            run_bsqr(row, projects_home, task)
        finally:
            fix_file_permissions(projects_home, row)


def run_bsqr(row, projects_home, task):
    log_info("Running Base quality score recalibration job (SAMPLE_ID='%s', BSQR='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.bsqr(), task, row.genomebuild()))

    if row.bsqr() == "no-bsqr":
        log_info("[%s] Skipping Base quality score recalibration for sample: '%s'" % (row.bsqr(), row.sample_id()))
        return

    bsqr_data = sample.bsqr_data(row, projects_home)

    if os.path.isfile(bsqr_data.bsqr_bam_out()):
        log_info("Already run bam recab..Skipping %s" % bsqr_data.bsqr_bam_out())
        return

    if not os.path.isfile(bsqr_data.bsqr_bam_in()):
        raise IOError("Can not find required BAM file: %s" % bsqr_data.bsqr_bam_in())

    genomebuild = select_genomebuild(row, projects_home)

    callables = {
        "bam-bsqr|no-task": bam_bsqr,
        "gatk-bsqr|no-task": gatk_bsqr,
    }

    callables.get(row.aligner() + "|" + task, unrecognized_options)(row, task, bsqr_data, genomebuild)


def unrecognized_options(row, task, *args):
    if callable is None:
        raise ValueError(
            "Unrecognised BSQR options (BSQR='%s', TASK='%s')" % (row.bsqr(), task))


def select_genomebuild(row, projects_home):
    genomebuild = genome_build.select(row.genomebuild(), projects_home)
    if genomebuild.known_indels() is None:
        raise ValueError(
            "No genome selected for sample '%s'. GENOMEBUILD value is '%s', but it should be one of [b37, hg19]" %
            row.genomebuild())
    return genomebuild


def bam_bsqr(row, task, bsqr_data, genomebuild):
    log_debug("bam_bsqr (SAMPLE_ID='%s', BSQR='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.bsqr(), task, row.genomebuild()))

    run_script("bam-bsqr", "bam-bsqr.tmpl.sh",
               NCPU=str(row.ncpu()),
               TMP_DIR=bsqr_data.tmp_dir(),
               BAM_IN=bsqr_data.bsqr_bam_in(),
               BAM_OUT=bsqr_data.bsqr_bam_out(),
               REFFASTA=genomebuild.ref_fasta(),
               DBSNP_RECAB=genomebuild.dbsnp_recab())


def gatk_bsqr(row, task, bsqr_data, genomebuild):
    log_debug("gatk_bsqr (SAMPLE_ID='%s', BSQR='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.bsqr(), task, row.genomebuild()))

    run_script("gatk-bsqr", "gatk-bsqr.tmpl.sh",
               NCPU=str(row.ncpu()),
               TMP_DIR=bsqr_data.tmp_dir(),
               BAM_IN=bsqr_data.bsqr_bam_in(),
               BAM_OUT=bsqr_data.bsqr_bam_out(),
               REFFASTA=genomebuild.ref_fasta(),
               DBSNP_RECAB=genomebuild.dbsnp_recab(),

               KNOWN_INDELS=genomebuild.known_indels(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               KNOWN_SNPS_OMNI=genomebuild.known_snps_omni(),
               KNOWN_SNPS_1000G=genomebuild.known_snps_1000g(),
               RECAL_DATA_TABLE=bsqr_data.reports_path(bsqr_data.bam_prefix() + ".recal_data.table"))


def run_script(dir, filename, **kwargs):
    tmpl = sh_template.load("bsqr", dir, filename)
    run_command(tmpl.create_sh_file(**kwargs))


if __name__ == '__main__':
    main(sys.argv[1:])
