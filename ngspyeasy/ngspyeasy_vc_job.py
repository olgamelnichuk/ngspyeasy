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
from logger import init_logger, logger

LOGGER_NAME = "variant_calling"

BASE_QUAL = "20"
MAP_QUAL = "20"
COVERAGE_MIN = "2"


def fix_file_permissions(projects_home, row, uid, gid):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id(), uid, gid)


def main(argv):
    args = cmdargs.parse_job_args(argv, "Variant Calling")

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    init_logger(log_file, args.verbose)
    logger().info("Starting...")
    logger().debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    logger().info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        logger().error(e)
        sys.exit(1)

    try:
        ngspyeasy_vc_job(tsv_conf, projects_home, args.sample_id, args.task, args.uid, args.gid)
    except Exception as ex:
        logger().exception(ex)
        sys.exit(1)


def ngspyeasy_vc_job(tsv_conf, projects_home, sample_id, task, uid, gid):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        try:
            run_vc(row, projects_home, task)
        finally:
            fix_file_permissions(projects_home, row, uid, gid)


def run_vc(row, projects_home, task):
    logger().info("Running variant caller job (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    if row.varcaller() == "no-vc":
        logger().info("[%s] Skipping variant calling for sample: '%s'" % (row.varcaller(), row.sample_id()))
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

    logger().info("Genome build selected: '%s'" % genomebuild.refdir())
    return genomebuild


def filter_duplicates(ngs_type):
    return 1 if ngs_type == "TGS" or ngs_type == "WEX" else 0


def prepare(row, projects_home, task):
    logger().debug("prepare (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    bam_file = vc_data.vc_bam_in()

    if not os.path.isfile(bam_file):
        raise IOError("Can not find [%s] for Variant Calling." % bam_file)

    mapped_reads_bed = vc_data.reports_path(os.path.basename(bam_file) + ".mapped.reads.bed")
    genomecov_bed = vc_data.reports_path(os.path.basename(bam_file) + ".genomecov.bed")
    logger().info("Mapped Reads: %s" % mapped_reads_bed)
    logger().info("Callable Regions: %s" % genomecov_bed)

    if not os.path.isfile(mapped_reads_bed) or not os.path.isfile(genomecov_bed):
        logger().info("Mapped Reads BED File and Callable Regions File do not exist. Generating...")
        run_script(tmpl_path("prepare", "callable-regions"),
                   NCPU=row.ncpu(),
                   BAM_FILE=bam_file,
                   MAPPED_READS_BED=mapped_reads_bed,
                   GENOMECOV_BED=genomecov_bed)

    filtered_bam_file = vc_data.vc_filtered_bam()
    logger().info("Filtered BAM file: %s" % filtered_bam_file)

    if not os.path.isfile(filtered_bam_file):
        logger().info("Filtered BAM file doesn't exist. Filtering BAM... (Q20 and proper_pair)")
        run_script(tmpl_path("prepare", "filter-bam"),
                   NCPU=row.ncpu(),
                   BAM_FILE=bam_file,
                   FILTERED_BAM=filtered_bam_file)


def freebayes_parallel(row, projects_home, task):
    logger().debug("freebayers_parallel (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using freebayes-parallel..." % vcf_gz)
    run_script(tmpl_path(row.varcaller()),
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
    logger().debug("freebayers_default (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using freebayes-default..." % vcf_gz)
    run_script(tmpl_path(row.varcaller()),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz())


def platypus(row, projects_home, task):
    logger().debug("platypus (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using platypus..." % vcf_gz)
    run_script(tmpl_path(row.varcaller()),
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
    logger().debug("platypus-default (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using platypus-default..." % vcf_gz)
    run_script(tmpl_path(row.varcaller()),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               FILTER_DUPLICATES=filter_duplicates(row.ngs_type()),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz())


def unified_genotyper(row, projects_home, task):
    logger().debug("unified-genotyper (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using unified-genotyper..." % vcf_gz)
    run_script(tmpl_path(row.varcaller()),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz(),
               TMP_DIR=vc_data.tmp_dir())


def haplotype_caller(row, projects_home, task):
    logger().debug("haplotype_caller (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using haplotype-caller..." % vcf_gz)
    run_script(tmpl_path(row.varcaller()),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               RAW_VCF=vc_data.raw_vcf(),
               VCF_GZ=vc_data.vcf_gz(),
               RAW_VCF_GZ=vc_data.raw_vcf_gz(),
               TMP_DIR=vc_data.tmp_dir())


def ensemble_freebayes_parallel(row, projects_home, task):
    logger().debug("ensemble_freebayers_parallel (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz("freebayers_parallel")

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using ensemble.freebayes-parallel..." % vcf_gz)
    run_script(tmpl_path(row.varcaller(), "freebayers"),
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
    logger().debug("ensemble_platypus (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz("platypus")

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using ensemble.platypus..." % vcf_gz)
    run_script(tmpl_path(row.varcaller(), "platypus"),
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
    logger().debug("ensemble_haplotype_caller (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz("HaplotypeCaller")

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using haplotype-caller..." % vcf_gz)
    run_script(tmpl_path(row.varcaller(), "HaplotypeCaller"),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               FILTERED_BAM=vc_data.vc_filtered_bam(),
               KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
               RAW_VCF=vc_data.raw_vcf("HaplotypeCaller"),
               VCF_GZ=vc_data.vcf_gz("HaplotypeCaller"),
               RAW_VCF_GZ=vc_data.raw_vcf_gz("HaplotypeCaller"),
               TMP_DIR=vc_data.tmp_dir())


def ensemble_bcbio_variation(row, projects_home, task):
    logger().debug("ensemble_bcbio_variation (SAMPLE_ID='%s', VARCALLER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.varcaller(), task, row.genomebuild()))

    vc_data = sample.vc_data(row, projects_home)
    vcf_gz = vc_data.vcf_gz()

    logger().info("VCF file: %s" % vcf_gz)

    if os.path.isfile(vcf_gz):
        logger().info("VCF file (%s) exists. Skipping this bit.." % vcf_gz)
        return

    genomebuild = select_genomebuild(row, projects_home)

    logger().info("VCF file (%s) doesn't exist. Running Variant Calling using bcbio-variation-recall..." % vcf_gz)
    run_script(tmpl_path(row.varcaller(), "bcbio-variation"),
               NCPU=row.ncpu(),
               REFFASTA=genomebuild.ref_fasta(),
               VCF_GZ=vc_data.vcf_gz(),
               FREEBAYERS_VCF_GZ=vc_data.vcf_gz("freebayes-parallel"),
               HAPLOTYPE_CALLER_VCF_GZ=vc_data.vcf_gz("HaplotypeCaller"))


def tmpl_path(dir, filename=None):
    return os.path.join(dir, (dir if filename is None else filename) + ".tmpl.sh")


def run_script(path, **kwargs):
    tmpl = sh_template.load("vc", path)
    run_command(tmpl.create_sh_file(**kwargs))


if __name__ == '__main__':
    main(sys.argv[1:])
