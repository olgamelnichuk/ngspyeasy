#!/usr/bin/env python
import datetime
from signal import signal, SIGPIPE, SIG_DFL
import subprocess
import sys

import cmdargs
from shcmd import run_command
import sh_template
import os
import sample
import genome_build
import projects_dir
import tsv_config
from logger import init_logger, logger


def fix_file_permissions(projects_home, row):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id())


def main(argv):
    args = cmdargs.parse_job_args(argv, "Alignment")

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
        ngspyeasy_alignment_job(tsv_conf, projects_home, args.sample_id, args.task)
    except Exception as ex:
        logger().exception(ex)
        sys.exit(1)


def ngspyeasy_alignment_job(tsv_conf, projects_home, sample_id, task=None):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = [x for x in rows2run if x.sample_id() == sample_id]

    for row in rows2run:
        try:
            run_alignment(row, projects_home, task)
        finally:
            fix_file_permissions(projects_home, row)


def run_alignment(row, projects_home, task):
    logger().info("Running alignmnent job (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    if row.aligner() == "no-align":
        logger().info("[%s] Skipping alignment" % row.aligner())
        return

    callables = {
        "bwa|no-task": bwa,
        "novoalign|no-task": novoalign,
        "stampy|bwa": stampy_bwa,
        "stampy|stampy": stampy_stampy,
        "stampy|picard_cleansam": stampy_cleansam,
        "stampy|picard_addorreplacereadgroups": stampy_picard_addorreplacereadgroups,
        "bowtie2|no-task": bowtie2,
        "snap|no-task": snap
    }

    callables.get(row.aligner() + "|" + task, unrecognized_options)(row, projects_home, task)


def unrecognized_options(row, projects_home, task):
    if callable is None:
        raise ValueError(
            "Unrecognised options (ALIGNER='%s', TASK='%s')" % (row.aligner(), task))


def alignment_results_exist(align_data):
    bam_out = align_data.dupl_mark_bam()
    bam_bed = align_data.dupl_mark_bed()
    exist = os.path.isfile(bam_out) and os.path.isfile(bam_bed)
    if exist:
        logger().info("Looks like aligned BAM file already exists [%s, %s] Skipping Alignment..." % (bam_out, bam_bed))
    return exist


def select_genome(row, projects_home):
    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome is None:
        raise ValueError("No genome selected. Choose one of [b37] or [hg19]")
    return genome


def common_parameters(align_data, projects_home):
    fastq_files = align_data.trimmed_fastq_files()
    not_exist = [x for x in fastq_files if not os.path.isfile(x)]
    if len(not_exist) > 0:
        raise ValueError("Required FastQ files do not exist: %s" % not_exist)

    row = align_data.row()
    return dict(
        NCPU=row.ncpu(),
        BAM_PREFIX=align_data.bam_prefix(),
        PLATFORM_UNIT=find_platform_unit(fastq_files[0]),
        NGS_PLATFORM=row.ngs_platform(),
        DNA_PREP_LIBRARY_ID=row.dna_prep_library_id(),
        RUNDATE=datetime.datetime.now().strftime("%d%m%y%H%M%S"),
        FQ1=fastq_files[0],
        FQ2=fastq_files[1],
        GENOMEINDEX=select_genome(row, projects_home).genome_index(),
        DISCORDANT_SAM=align_data.discordant_sam(),
        DISCORDANT_BAM=align_data.discordant_bam(),
        SPLITREAD_SAM=align_data.splitread_sam(),
        SPLITREAD_BAM=align_data.splitread_bam(),
        UNMAPPED_FASTQ=align_data.unmapped_fastq(),
        DUPEMARK_BAM=align_data.dupl_mark_bam(),
        DUPEMARK_FLAGSTAT_REPORT=align_data.dupl_mark_bam_flagstat(),
        DUPEMARK_BED_REPORT=align_data.dupl_mark_bed(),
        TMP_DIR=align_data.tmp_dir(),
        ALIGNMENTS_DIR=align_data.alignments_dir()
    )


def bwa(row, projects_home, task):
    logger().debug("bwa (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    run_script("bwa", "bwa.tmpl.sh", **common_parameters(align_data, projects_home))


def novoalign(row, projects_home, task):
    logger().debug("novoalign (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    params = common_parameters(align_data, projects_home)
    params.update(
        K_STATS=align_data.alignments_path(align_data.bam_prefix() + ".K.stats"))

    run_script("novoalign", "novoalign.tmpl.sh", **params)


def bowtie2(row, projects_home, task):
    logger().debug("bowtie2 (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    params = common_parameters(align_data, projects_home)
    params.update(
        FAST="--end-to-end --sensitive",
        SLOW="--local --sensitive-local")

    run_script("bowtie2", "bowtie2.tmpl.sh", **params)


def snap(row, projects_home, task):
    logger().debug("snap (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    params = common_parameters(align_data, projects_home)
    params.update(REFDIR=select_genome(row, projects_home).refdir())

    run_script("snap", "snap.tmpl.sh", **params)


def stampy_bwa(row, projects_home, task):
    logger().debug("stampy_bwa (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    fastq_files = align_data.trimmed_fastq_files()
    tmp_bam = align_data.tmp_bam()

    if os.path.isfile(tmp_bam):
        logger().info("tmp.bam already exists (%s)... skipping bwa task" % tmp_bam)
        return

    run_script("stampy", "bwa.tmpl.sh",
               NCPU=row.ncpu(),
               BAM_PREFIX=align_data.bam_prefix(),
               PLATFORM_UNIT=find_platform_unit(fastq_files[0]),
               NGS_PLATFORM=row.ngs_platform(),
               DNA_PREP_LIBRARY_ID=row.dna_prep_library_id(),
               RUNDATE=datetime.datetime.now().strftime("%d%m%y%H%M%S"),
               FQ1=fastq_files[0],
               FQ2=fastq_files[1],
               GENOMEINDEX=select_genome(row, projects_home).genome_index(),
               TMP_BAM=tmp_bam,
               TMP_DIR=align_data.tmp_dir(),
               ALIGNMENTS_DIR=align_data.alignments_dir())


def stampy_stampy(row, projects_home, task):
    logger().debug("stampy_stampy (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    if os.path.isfile(align_data.dupl_mark_tmp_bam()):
        logger().info("dupemk.tmp.bam already exists (%s)... skipping stampy task" % align_data.dupl_mark_tmp_bam())
        return

    if not os.path.isfile(align_data.tmp_bam()):
        raise IOError("Can't find BWA input for stampy: %s" % align_data.tmp_bam())

    run_script("stampy", "stampy.tmpl.sh",
               STAMPY_VERSION="stampy-1.0.27",
               NCPU=row.ncpu(),
               GENOMEINDEX=select_genome(row, projects_home).genome_index(),
               DISCORDANT_SAM=align_data.discordant_sam(),
               DISCORDANT_BAM=align_data.discordant_bam(),
               SPLITREAD_SAM=align_data.splitread_sam(),
               SPLITREAD_BAM=align_data.splitread_bam(),
               UNMAPPED_FASTQ=align_data.unmapped_fastq(),
               DUPEMARK_TMP_BAM=align_data.dupl_mark_tmp_bam(),
               TMP_BAM=align_data.tmp_bam(),
               TMP_DIR=align_data.tmp_dir(),
               ALIGNMENTS_DIR=align_data.alignments_dir())


def stampy_cleansam(row, projects_home, task):
    logger().debug("stampy_cleansam (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    if os.path.isfile(align_data.dupl_mark_tmp_cleansam_bam()):
        logger().info(
            "dupemk.tmpcleansam.bam already exists (%s)... skipping picard-cleansam task" % align_data.dupl_mark_tmp_cleansam_bam())
        return

    if not os.path.isfile(align_data.dupl_mark_tmp_bam()):
        raise IOError("Can't find input for picard-cleansam: %s" % align_data.dupl_mark_tmp_bam())

    run_script("stampy", "picard_cleansam.tmpl.sh",
               PICARD_VERSION="picard-tools-1.128",
               DUPEMARK_TMP_BAM=align_data.dupl_mark_tmp_bam(),
               DUPEMARK_CLEANSAM_BAM=align_data.dupl_mark_tmp_cleansam_bam(),
               TMP_DIR=align_data.tmp_dir())


def stampy_picard_addorreplacereadgroups(row, projects_home, task):
    logger().debug(
        "stampy_picard_addorreplacereadgroups (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    align_data = sample.alignment_data(row, projects_home)
    if alignment_results_exist(align_data):
        return

    if not os.path.isfile(align_data.dupl_mark_tmp_cleansam_bam()):
        raise IOError(
            "Can't find input for picard-addorreplacereadgroups task" % align_data.dupl_mark_tmp_cleansam_bam())

    fastq_files = align_data.trimmed_fastq_files()
    run_script("stampy", "picard_addorreplacereadgroups.tmpl.sh",
               PICARD_VERSION="picard-tools-1.128",
               NCPU=row.ncpu(),
               BAM_PREFIX=align_data.bam_prefix(),
               PLATFORM_UNIT=find_platform_unit(fastq_files[0]),
               NGS_PLATFORM=row.ngs_platform(),
               DNA_PREP_LIBRARY_ID=row.dna_prep_library_id(),
               RUNDATE=datetime.datetime.now().strftime("%d%m%y%H%M%S"),
               DUPEMARK_BAM=align_data.dupl_mark_bam(),
               DUPEMARK_FLAGSTAT_REPORT=align_data.dupl_mark_bam_flagstat(),
               DUPEMARK_BED_REPORT=align_data.dupl_mark_bed(),
               DUPEMARK_CLEANSAM_BAM=align_data.dupl_mark_tmp_cleansam_bam(),
               DUPEMARK_TMP_BAM=align_data.dupl_mark_tmp_bam(),
               TMP_DIR=align_data.tmp_dir())


def run_script(dir, filename, **kwargs):
    tmpl = sh_template.load("alignment", dir, filename)
    run_command(tmpl.create_sh_file(**kwargs))


def find_platform_unit(fastq_file):
    cmd = "zcat %s | head -1 | sed 's/:/\\t/' - | cut -f 1 | sed 's/@//g' - " % fastq_file
    logger().debug("platform_info=[%s]" % cmd)
    out = subprocess.check_output(cmd, shell=True, preexec_fn=lambda: signal(SIGPIPE, SIG_DFL))
    return out.strip()


if __name__ == '__main__':
    main(sys.argv[1:])
