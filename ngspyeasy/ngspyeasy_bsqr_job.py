#!/usr/bin/env python
import sys

from shutils import script_from_template, run_command
import os
import cmdargs
import projects_dir
import tsv_config
import sample_data
import genome_build
from logger import init_logger, get_logger

LOGGER_NAME = "bsqr"


def log_info(msg):
    get_logger(LOGGER_NAME).info(msg)


def log_debug(msg):
    get_logger(LOGGER_NAME).debug(msg)


def log_error(msg):
    get_logger(LOGGER_NAME).error(msg)


def main(argv):
    args = cmdargs.parse_job_args(argv, "Base quality score recalibration")

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
        ngspyeasy_bsqr_job(tsv_conf, projects_home, args.sample_id)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_bsqr_job(tsv_conf, projects_home, sample_id):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_bsqr(row, projects_home)


def run_bsqr(row, projects_home):
    log_info("Running Base quality score recalibration job (SAMPLE_ID='%s', BSQR='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.bsqr(), row.genomebuild()))

    if row.bsqr() not in ["gatk-bsqr", "bam-bsqr", "no-bsqr"]:
        raise ValueError(
            "Unrecognised BSQR option. Should be one of [bam-bsqr] [gatk-bsqr] or [no-bsqr]: '%s'" % row.realn())

    if row.bsqr() == "no-bsqr":
        log_info("[%s] Skipping Base quality score recalibration for sample: '%s'" % (row.bsqr(), row.sample_id()))
        return

    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome.known_indels() is None:
        raise ValueError(
            "No genome selected for sample '%s'. GENOMEBUILD value is '%s', but it should be one of [b37, hg19]" %
            row.genomebuild())

    log_info("Genome build selected: '%s'" % genome.refdir())

    sample = sample_data.create(row, projects_home)

    base_dir = os.path.dirname(__file__)
    template_path = os.path.join(base_dir, "resources", "bsqr", row.bsqr(), row.bsqr() + ".tmpl.sh")

    log_debug("Using script template file: %s" % template_path)

    script = script_from_template(template_path)

    log_debug("Script template to run: %s" % script.source())

    script.add_variables(
        NCPU=str(row.ncpu()),
        REFFASTA=genome.ref_fasta(),
        DBSNP_RECAB=genome.dbsnp_recab(),
        TMP_DIR=sample.tmp_dir()
    )

    if row.bsqr() == "bam-bsqr":
        if os.path.isfile(sample.dupl_mark_realn_bam("bam-realn")):
            bam_in = sample.dupl_mark_realn_bam("bam-realn")
            bam_out = sample.dupl_mark_realn_bsqr_bam("bam-realn")
        elif os.path.isfile(sample.dupl_mark_bam()):
            bam_in = sample.dupl_mark_bam()
            bam_out = sample.dupl_mark_bsqr_bam()
        else:
            raise IOError("Can not find required BAM files in %s" % sample.alignments_dir())

        log_info("BAM in: %s" % bam_in)
        log_info("BAM out: %s" % bam_out)

        if os.path.isfile(bam_out):
            log_info("Already run bam recab..Skipping")
            return

        script.add_variables(
            BAM_IN=bam_in,
            BAM_OUT=bam_out
        )

    elif row.bsqr() == "gatk-bsqr":
        if os.path.isfile(sample.dupl_mark_realn_bam("gatk-realn")):
            bam_in = sample.dupl_mark_realn_bam("gatk-realn")
            bam_out = sample.dupl_mark_realn_bsqr_bam("gatk-realn")
        elif os.path.isfile(sample.dupl_mark_bam()):
            bam_in = sample.dupl_mark_bam()
            bam_out = sample.dupl_mark_bsqr_bam()
        else:
            raise IOError("Can not find required BAM files in %s" % sample.alignments_dir())

        log_info("BAM in: %s" % bam_in)
        log_info("BAM out: %s" % bam_out)

        if os.path.isfile(bam_out):
            log_info("Already run GATK BSQR..Skipping")
            return

        script.add_variables(
            BAM_IN=bam_in,
            BAM_OUT=bam_out,
            KNOWN_INDELS=genome.known_indels(),
            KNOWN_SNPS_b138=genome.known_snps_b138(),
            KNOWN_SNPS_OMNI=genome.known_snps_omni(),
            KNOWN_SNPS_1000G=genome.known_snps_1000g(),
            RECAL_DATA_TABLE=sample.reports_path(sample.bam_prefix() + ".recal_data.table")
        )

    log_debug("Script template variables:\n %s" % "\n".join(script.variable_assignments()))

    run_command(script.to_temporary_file(), get_logger(LOGGER_NAME))

if __name__ == '__main__':
    main(sys.argv[1:])
