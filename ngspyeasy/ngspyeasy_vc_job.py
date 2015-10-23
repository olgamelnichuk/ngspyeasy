#!/usr/bin/env python
import sys

import sh_template

from shcmd import run_command

import os
import cmdargs
import projects_dir
import tsv_config
import genome_build
import sample
from logger import init_logger, log_info, log_debug, log_error, log_exception

LOGGER_NAME = "variant_calling"

BASE_QUAL = "20"
MAP_QUAL = "20"
COVERAGE_MIN = "2"


def fix_file_permissions(projects_home, row):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id())


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
        log_exception(ex)
        sys.exit(1)


def ngspyeasy_vc_job(tsv_conf, projects_home, sample_id, task):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        try:
            run_vc(row, projects_home, task)
        finally:
            fix_file_permissions(projects_home, row)


def run_vc(row, projects_home, task):
    log_info("Running variant caller job (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    if row.varcaller() == "no-vc":
        log_info("[%s] Skipping variant calling for sample: '%s'" % (row.varcaller(), row.sample_id()))
        return

    callables = {
        "freebayes-parallel|no-task": freebayes_parallel,
        "freebayes-default|no-task": freebayes_default,
        "platypus|no-task": platypus,
        "platypus-default|no-task": platypus_default,
        "UnifiedGenotyper|no-task": unified_genotyper,
        "HaplotypeCaller|no-task": haplotype_caller,
        "ensemble|freebayes-parallel": ensemble_freebayes_parallel,
        "ensemble|platypus": ensemble_platypus,
        "ensemble|HaplotypeCaller": ensemble_haplotype_caller,
        "ensemble|bcbio-variation": ensemble_bcbio_variation
    }

    if task == "prepare":
        prepare(row, projects_home, task)
    else:
        callables.get(row.varcaller() + "|" + task, unrecognized_options)(row, projects_home, task)


def unrecognized_options(row, projects_home, task):
    if callable is None:
        raise ValueError(
            "Unrecognised options (VARCALLER='%s', TASK='%s')" % (row.varcaller(), task))


def select_genomebuild(row, projects_home):
    genomebuild = genome_build.select(row.genomebuild(), projects_home)
    if genomebuild.known_indels() is None:
        raise ValueError(
            "No genome selected for sample '%s'. GENOMEBUILD value is '%s', but it should be one of [b37, hg19]" % (
                row.sample_id(), row.genomebuild()))

    log_info("Genome build selected: '%s'" % genomebuild.refdir())
    return genomebuild


def filter_duplicates(ngs_type):
    return 1 if ngs_type == "TGS" or ngs_type == "WEX" else 0


def prepare(row, projects_home, task):
    log_debug("prepare (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    bam_file = vc_data.vc_bam_in()

    if not os.path.isfile(bam_file):
        raise IOError("Can not find [%s] for Variant Calling." % bam_file)

    mapped_reads_bed = vc_data.reports_path(os.path.basename(bam_file) + ".mapped.reads.bed")
    genomecov_bed = vc_data.reports_path(os.path.basename(bam_file) + ".genomecov.bed")
    log_info("Mapped Reads: %s" % mapped_reads_bed)
    log_info("Callable Regions: %s" % genomecov_bed)

    if not os.path.isfile(mapped_reads_bed) or not os.path.isfile(genomecov_bed):
        log_info("Mapped Reads BED File and Callable Regions File do not exist. Generating...")
        run_script("prepare", "callable-regions.tmpl.sh",
                   NCPU=row.ncpu(),
                   BAM_FILE=bam_file,
                   MAPPED_READS_BED=mapped_reads_bed,
                   GENOMECOV_BED=genomecov_bed)

    filtered_bam_file = vc_data.vc_filtered_bam()
    log_info("Filtered BAM file: %s" % filtered_bam_file)

    if not os.path.isfile(filtered_bam_file):
        log_info("Filtered BAM file doesn't exist. Filtering BAM... (Q20 and proper_pair)")
        run_script("prepare", "filter-bam.tmpl.sh",
                   NCPU=row.ncpu(),
                   BAM_FILE=bam_file,
                   FILTERED_BAM=filtered_bam_file)


def freebayes_parallel(row, projects_home, task):
    log_debug("freebayers_parallel (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using freebayes-parallel..." % vcf_gz)
    run_script(row.varcaller(), row.varcaller(),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               COVERAGE_MIN=COVERAGE_MIN,
               MAP_QUAL=MAP_QUAL,
               BASE_QUAL=BASE_QUAL,
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz())


def freebayes_default(row, projects_home, task):
    log_debug("freebayers_default (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using freebayes-default..." % vcf_gz)
    run_script(row.varcaller(), row.varcaller(),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz())


def platypus(row, projects_home, task):
    log_debug("platypus (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using platypus..." % vcf_gz)
    run_script(row.varcaller(), row.varcaller(),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               FILTER_DUPLICATES=filter_duplicates(row.ngs_type()),
               COVERAGE_MIN=COVERAGE_MIN,
               MAP_QUAL=MAP_QUAL,
               BASE_QUAL=BASE_QUAL,
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz())


def platypus_default(row, projects_home, task):
    log_debug("platypus-default (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using platypus-default..." % vcf_gz)
    run_script(row.varcaller(), row.varcaller(),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               FILTER_DUPLICATES=filter_duplicates(row.ngs_type()),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz())


def unified_genotyper(row, projects_home, task):
    log_debug("unified-genotyper (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using unified-genotyper..." % vcf_gz)
    run_script(row.varcaller(), row.varcaller(),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz(),
               TMP_DIR=vc_data.tmp_dir())


def haplotype_caller(row, projects_home, task):
    log_debug("haplotype_caller (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using haplotype-caller..." % vcf_gz)
    run_script(row.varcaller(), row.varcaller(),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz(),
               TMP_DIR=vc_data.tmp_dir())


def ensemble_freebayes_parallel(row, projects_home, task):
    log_debug("ensemble_freebayers_parallel (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz("freebayers_parallel")

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using ensemble.freebayes-parallel..." % vcf_gz)
    run_script(row.varcaller(), "freebayers",
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               COVERAGE_MIN=COVERAGE_MIN,
               MAP_QUAL=MAP_QUAL,
               BASE_QUAL=BASE_QUAL,
               RAW_VCF=vc_data.raw_vcf("freebayers_parallel"),
               VCF_GZ=vc_data.vcf_gz("freebayers_parallel"),
               RAW_VCF_GZ=vc_data.raw_vcf_gz("freebayers_parallel"))


def ensemble_platypus(row, projects_home, task):
    log_debug("ensemble_platypus (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz("platypus")

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using ensemble.platypus..." % vcf_gz)
    run_script(row.varcaller(), "platypus",
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               FILTER_DUPLICATES=filter_duplicates(row.ngs_type()),
               COVERAGE_MIN=COVERAGE_MIN,
               MAP_QUAL=MAP_QUAL,
               BASE_QUAL=BASE_QUAL,
               RAW_VCF=vc_data.raw_vcf("platypus"),
               VCF_GZ=vc_data.vcf_gz("platypus"),
               RAW_VCF_GZ=vc_data.raw_vcf_gz("platypus"))


def ensemble_haplotype_caller(row, projects_home, task):
    log_debug("ensemble_haplotype_caller (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz("HaplotypeCaller")

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using haplotype-caller..." % vcf_gz)
    run_script(row.varcaller(), "HaplotypeCaller",
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               RAW_VCF=vc_data.raw_vcf("HaplotypeCaller"),
               VCF_GZ=vc_data.vcf_gz("HaplotypeCaller"),
               RAW_VCF_GZ=vc_data.raw_vcf_gz("HaplotypeCaller"),
               TMP_DIR=vc_data.tmp_dir())


def ensemble_bcbio_variation(row, projects_home, task):
    log_debug("ensemble_bcbio_variation (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    log_info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        log_info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    log_info("VCF file (%s) doesn't exist. Running Variant Calling using bcbio-variation-recall..." % vcf_gz)
    run_script(row.varcaller(), "bcbio-variation",
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               VCF_GZ=vc_data.vcf_gz(),
               FREEBAYERS_VCF_GZ=vc_data.vcf_gz("freebayes-parallel"),
               HAPLOTYPE_CALLER_VCF_GZ=vc_data.vcf_gz("HaplotypeCaller"))


def run_script(dir, filename, **kwargs):
    tmpl = sh_template.load("vc", dir, filename)
    run_command(tmpl.create_sh_file(**kwargs))


if __name__ == '__main__':
    main(sys.argv[1:])
