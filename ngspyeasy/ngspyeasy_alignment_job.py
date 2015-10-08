#!/usr/bin/env python
import getopt
import datetime
from signal import signal, SIGPIPE, SIG_DFL
import subprocess
import sys

from shutils import script_from_template, run_command
import os
import sample_data
import projects_dir
import tsv_config
from cmdline_options import check_cmdline_options
from logger import init_logger, log_error, log_set_current_step, log_info, log_debug


def usage():
    print """
Usage:  ngspyeasy_alignment_job -c <config_file> -d <project_directory> -i <sample_id> -t <task>

Options:
        -c  STRING  configuration file
        -d  STRING  project directory
        -v  NULL    verbose
        -h  NULL    show this message
        -i  STRING sample id
        -t  STRING task name (e.g. stampy aligner tasks are: ['bwa', 'stampy', 'picard_cleansam', 'picard_addorreplacereadgroups'])
"""


def exit_with_error(msg):
    print >> sys.stderr, "ERROR:" + msg
    sys.exit(1)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hvc:d:i:t:", ["help"])
        if len(opts) == 0:
            usage()
            sys.exit(1)

    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    tsv_config_file = None
    ngs_projects_dir = None
    verbose = False
    sample_id = None
    task = None
    for opt, val in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == "-c":
            tsv_config_file = val
        elif opt == "-d":
            ngs_projects_dir = val
        elif opt == "-v":
            verbose = True
        elif opt == "-i":
            sample_id = val
        elif opt == "-t":
            task = val
        else:
            assert False, "unhandled option"

    (tsv_name, projects_home, errmsg) = check_cmdline_options(tsv_config_file, ngs_projects_dir)
    if errmsg:
        exit_with_error(errmsg)

    init_logger(projects_dir.sample_log_file(projects_home, tsv_name, sample_id), verbose)
    log_set_current_step("ngspyeasy_fastqc_job")

    tsv_conf = tsv_config.parse(projects_dir.config_full_path(projects_home, tsv_name))
    if tsv_conf is None:
        exit_with_error("Invalid TSV config. See logs for details...")

    try:
        ngspyeasy_alignment_job(tsv_conf, projects_home, sample_id, task)
    except Exception as ex:
        log_error(ex)
        sys.exit(1)


def ngspyeasy_alignment_job(tsv_conf, projects_home, sample_id, task=None):
    rows2run = tsv_conf.all_rows()
    if sample_id is not None:
        rows2run = filter(lambda x: x.sample_id() == sample_id, rows2run)

    for row in rows2run:
        run_alignment(row, projects_home, task)


def run_alignment(row, projects_home, task):
    if task is None:
        task = row.aligner()

    log_info("Running alignmnent job (SAMPLE_ID='%s', ALIGNER='%s', TASK='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), task, row.genomebuild()))

    if row.aligner() == "no-align":
        log_info("[%s] Skipping alignment" % row.aligner())
        return

    if row.aligner() not in ["bwa", "novoalign", "stampy", "bowtie2", "snap"]:
        raise ValueError("Unknown aligner value: '%s'" % row.aligner())

    ngs_resources = projects_dir.resources_dir(projects_dir)

    if row.genomebuild() == "b37":
        refdir = os.path.join(ngs_resources, "reference_genomes_b37")
        genome_index = os.path.join(refdir, "human_g1k_v37")
    elif row.genomebuild() == "hg19":
        refdir = os.path.join(ngs_resources, "reference_genomes_hg19")
        genome_index = os.path.join(refdir, "ucsc.hg19")
    elif row.genomebuild() == "hs37d5":
        refdir = os.path.join(ngs_resources, "reference_genomes_hs37d5")
        genome_index = os.path.join(refdir, "hs37d5")
    elif row.genomebuild() == "hs38DH":
        refdir = os.path.join(ngs_resources, "/reference_genomes_hs38DH")
        genome_index = os.path.join(refdir, "hs38DH")
    else:
        raise ValueError("No genome selected. Exiting. Choose one of [b37] or [hg19]")

    log_info("Genome build index: %s" % genome_index)
    log_info("Genome build dir: %s " % refdir)

    sample = sample_data.create(row, projects_home)

    if row.trim() == "no-trim":
        fastq = sample.fastq_files()
        log_info("TRIM set to '%s'. Using raw FastQ data: %s" % (row.trim(), fastq))
    elif row.trim in ["atrim", "btrim"]:
        fastq = sample.trimmomatic_paired_results()
        not_exist = [x for x in fastq if not os.path.isfile(x)]
        if len(not_exist) > 0:
            raise ValueError("Trimmed FastQC Data does not exist: %s" % not_exist)
        log_info("TRIM set to '%s'. Using trimmed FastQ data: %s" % (row.trim(), fastq))
    else:
        raise ValueError("Unknown TRIM value: '%s'" % row.trim())

    base_dir = os.path.dirname(__file__)
    template_path = os.path.join(base_dir, "resources", "alignment", row.aligner(), task + ".tmpl.sh")

    log_info("Using template file: %s" % template_path)

    script = script_from_template(template_path)

    log_info("Script template to run: %s" % script.source())

    bam_prefix = sample.bam_prefix()
    platform_unit = find_platform_unit(fastq[0])

    script.add_variables(
        NCPU=str(row.ncpu()),
        BAM_PREFIX=bam_prefix,
        PLATFORM_UNIT=platform_unit,
        NGS_PLATFORM=row.ngs_platform(),
        DNA_PREP_LIBRARY_ID=row.dna_prep_library_id(),
        RUNDATE=datetime.datetime.now().strftime("%d%m%y%H%M%S"),
        FQ1=fastq[0],
        FQ2=fastq[1],
        GENOMEINDEX=genome_index,
        DISCORDANT_SAM=sample.alignments_path(bam_prefix + ".discordant.sam"),
        DISCORDANT_BAM=sample.alignments_path(bam_prefix + ".discordant.bam"),
        SPLITREAD_SAM=sample.alignments_path(bam_prefix + ".splitread.sam"),
        SPLITREAD_BAM=sample.alignments_path(bam_prefix + ".splitread.bam"),
        UNMAPPED_FASTQ=sample.alignments_path(bam_prefix + ".unmapped.fastq"),
        DUPEMK_BAM=sample.alignments_path(bam_prefix + ".dupemk.bam"),
        DUPEMK_FLAGSTAT_REPORT=sample.reports_path(bam_prefix + ".dupemk.bam.flagstat"),
        DUPEMK_BED_REPORT=sample.reports_path(bam_prefix + ".dupemk.bed"),
        TMP_DIR=sample.tmp_dir(),
        ALIGNNMENTS_DIR=sample.alignments_dir()
    )

    if row.aligner() == "novoalign":
        script.add_variables(K_STATS=sample.alignments_path(bam_prefix + ".K.stats"))

    if row.aligner() == "stampy":
        script.add_variables(
            TMP_BAM=sample.alignments_path(bam_prefix + ".tmp.bam"),
            TMP_BAM_BAI=sample.alignments_path(bam_prefix + ".tmp.bam.bai"),
            DUPEMARK_TMP_BAM=sample.alignments_path(bam_prefix + ".dupemk.tmp.bam"),
            DUPEMARK_TMP_BAM_BAI=sample.alignments_path(bam_prefix + ".dupemk.tmp.bam.bai"),
            DUPEMARK_CLEANSAM_BAM=sample.alignments_path(bam_prefix + ".dupemk.tmpcleansam.bam"),
            DUPEMARK_CLEANSAM_BAM_BAI=sample.alignments_path(bam_prefix + ".dupemk.tmpcleansam.bam.bai"),
            STAMPY_VERSION="stampy-1.0.27",
            PICARD_VERSION="picard-tools-1.128"
        )

        if task not in ["bwa", "stampy", "picard_cleansam", "picard_addorreplacereadgroups"]:
            raise ValueError("Unknown stampy task: '%s'" % task)

        if task == "stampy":
            if not os.path.isfile(d["TMP_BAM"]):
                raise IOError("Tmp BAM file not found: %s" % d["TMP_BAM"])

    if row.aligner() == "bowtie2":
        script.add_variables(
            FAST="--end-to-end --sensitive",
            SLOW="--local --sensitive-local"
        )

    if row.aligner() == "snap":
        script.add_variables(REFDIR=refdir)

    log_debug("Script template variables:\n %s" % "\n".join(script.variable_assignments()))

    run_command(script.to_temporary_file(), logger)


def find_platform_unit(fastq_file):
    cmd = "zcat %s | head -1 | sed 's/:/\\t/' - | cut -f 1 | sed 's/@//g' - " % fastq_file
    log_debug("platform_info=[%s]" % cmd)
    out = subprocess.check_output(cmd, shell=True, preexec_fn=lambda: signal(SIGPIPE, SIG_DFL))
    return out.strip()


if __name__ == '__main__':
    main(sys.argv[1:])
