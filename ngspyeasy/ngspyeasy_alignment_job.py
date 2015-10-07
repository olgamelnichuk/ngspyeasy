#!/usr/bin/env python
import getopt
import datetime
from ngspyeasy import docker
import os
from signal import signal, SIGPIPE, SIG_DFL
import subprocess
import sys
import sample_data

import projects_dir
import tsv_config
from cmdline_options import check_cmdline_options
from logger import init_logger, log_error, log_set_current_step, log_info, log_debug

from string import Template


def usage():
    print """
Usage:  ngspyeasy_alignment_job -c <config_file> -d <project_directory> -i <sample_id> -t <task>

Options:
        -c  STRING  configuration file
        -d  STRING  project directory
        -v  NULL    verbose
        -h  NULL    show this message
        -i  STRING sample id
        -t  STRING task name (e.g. for stampy aligner tasks are: ['stampy_bwa', 'stampy_stampy', 'stampy_picard1', 'stampy_picard2'])
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
    log_info("Running alignmnent job (SAMPLE_ID='%s', ALIGNER='%s', GENOMEBUILD='%s')" % (
        row.sample_id(), row.aligner(), row.genomebuild()))

    if row.aligner() == "no-align":
        log_info("[%s] Skipping alignment" % row.aligner())
        return

    if row.genomebuild() == "b37":
        refdir = os.path.join(docker.NGS_RESOURCES, "reference_genomes_b37")
        genome_index = os.path.join(refdir, "human_g1k_v37")
    elif row.genomebuild() == "hg19":
        refdir = os.path.join(docker.NGS_RESOURCES, "reference_genomes_hg19")
        genome_index = os.path.join(refdir, "ucsc.hg19")
        ## hs37d5 and hs38DH added 04.08.15 generated using bwakit
    elif row.genomebuild() == "hs37d5":
        refdir = os.path.join(docker.NGS_RESOURCES, "reference_genomes_hs37d5")
        genome_index = os.path.join(refdir, "hs37d5")
    elif row.genomebuild() == "hs38DH":
        refdir = os.path.join(docker.NGS_RESOURCES, "/reference_genomes_hs38DH")
        genome_index = os.path.join(refdir, "hs38DH")
    else:
        raise ValueError("ERROR: No genome selected. Exiting. Choose one of [b37] or [hg19]")

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
            raise ValueError("ERROR:Trimmed FastQC Data does not exst: %s" % not_exist)
        log_info("TRIM set to '%s'. Using trimmed FastQ data: %s" % (row.trim(), fastq))

    bam_prefix = sample.bam_prefix()
    log_info("BAM prefix: '%s'" % bam_prefix)

    platform_unit = find_platform_unit(sample.fastq_files()[0])
    log_info("Platform unit: '%s'" % platform_unit)

    run_date = datetime.datetime.now().strftime("%d%m%y%H%M%S")

    discordant_sam = sample.alignments_path(bam_prefix + ".discordant.sam")
    discordant_bam = sample.alignments_path(bam_prefix + ".discordant.bam")
    splitter_sam = sample.alignments_path(bam_prefix + ".splitread.sam")
    splitter_bam = sample.alignments_path(bam_prefix + ".splitread.bam")
    unmapped_fastq = sample.alignments_path(bam_prefix + ".unmapped.fastq")
    dupemk_bam = sample.alignments_path(bam_prefix + ".dupemk.bam")
    dupemk_flagstat = sample.reports_path(bam_prefix + ".dupemk.bam.flagstat")
    dupemk_bed = sample.reports_path(bam_prefix + ".dupemk.bed")

    if row.aligner() == "bwa":

        shell_script = Template("""
time /usr/local/bin/bwa mem \
-M \
-t $NCPU \
-R '@RG\\tID:$BAM_PREFIX\tSM:$BAM_PREFIX\\tPU:$PLATFORM_UNIT}\\tPL:$NGS_PLATFORM\\tLB:$DNA_PREP_LIBRARY_ID\\tDT:$RUNDATE' \
$GENOMEINDEX.fasta \
$FQ1 $FQ2 | \
samblaster --addMateTags --excludeDups \
--discordantFile $DISCORDANT_SAM \
--splitterFile $SPLITTER_SAM \
--unmappedFile $UNMAPPED_FASTQ | \
sambamba view -t $NCPU -S -f bam /dev/stdin | \
sambamba sort -t $NCPU -m 2GB --tmpdir=$TMP_DIR -o $DUPEMK_BAM /dev/stdin && \
sambamba index $DUPEMK_BAM && \
sambamba flagstat -t $NCPU $DUPEMK_BAM > $DUPEMK_FLAGSTAT && \
bedtools bamtobed -i $DUPEMK_BAM | bedtools merge > $DUPEMK_BED && \
sambamba view -t $NCPU -S -f bam $DISCORDANT_SAM | \
sambamba sort -t $NCPU -m 2GB --tmpdir=$TMP_DIR -o $DISCORDANT_BAM /dev/stdin && \
sambamba index $DISCORDANT_BAM && \
sambamba view -t $NCPU -S -f bam $SPLITTER_SAM | \
sambamba sort -t $NCPU -m 2GB --tmpdir=$TMP_DIR -o $SPLITTER_BAM /dev/stdin && \
sambamba index $SPLITTER_BAM && \
rm -v $DISCORDANT_FILE && \
rm -v $SPLITTER_FILE && \
rm -rf $TMP_DIR/* && \
chmod -R 777 $ALIGNMENTS_DIR/*
        """)
        values = dict(NCPU=str(row.ncpu()),
                      BAM_PREFIX=bam_prefix,
                      PLATFORM_UNIT=platform_unit,
                      NGS_PLATFORM=row.ngs_platform(),
                      DNA_PREP_LIBRARY_ID=row.dna_prep_library_id(),
                      RUNDATE=run_date,
                      FQ1=fastq[0],
                      FQ2=fastq[1],
                      GENOMEINDEX=genome_index,
                      DISCORDANT_SAM=discordant_sam,
                      DISCORDANT_BAM=discordant_bam,
                      SPLITTER_SAM=splitter_sam,
                      SPLITTER_BAM=splitter_bam,
                      UNMAPPED_FASTQ=unmapped_fastq,
                      DUPEMK_BAM=dupemk_bam,
                      DUPEMK_FLAGSTAT=dupemk_flagstat,
                      DUPEMK_BED=dupemk_bed,
                      TMP_DIR=row.tmp_dir()
                      )


        log_info("\n" + shell_script.substitute(values))



def find_platform_unit(fastq_file):
    cmd = "zcat %s | head -1 | sed 's/:/\\t/' - | cut -f 1 | sed 's/@//g' - " % fastq_file
    log_debug("platform_info=[%s]" % cmd)
    out = subprocess.check_output(cmd, shell=True, preexec_fn=lambda: signal(SIGPIPE, SIG_DFL))
    return out.strip()


if __name__ == '__main__':
    main(sys.argv[1:])
