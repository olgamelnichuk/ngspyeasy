#!/usr/bin/env python

###
# Copyright 2015, EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###
import datetime
import os
from signal import signal, SIGPIPE, SIG_DFL
import subprocess
from logger import logger
import genome_build
import sample
from utils import LazyDict


def as_dict(row, projects_home):
    fqc_data = sample.fastqc_data(row, projects_home)
    align_data = sample.alignment_data(row, projects_home)
    realn_data = sample.realn_data(row, projects_home)
    bsqr_data = sample.bsqr_data(row, projects_home)
    vc_data = sample.vc_data(row, projects_home)

    trimmed_fq = align_data.trimmed_fastq_files()
    genomebuild = select_genome(row, projects_home)

    dict = LazyDict(
        NCPU=row.ncpu(),
        NGS_PLATFORM=row.ngs_platform(),
        DNA_PREP_LIBRARY_ID=row.dna_prep_library_id(),
        RUNDATE=datetime.datetime.now().strftime("%d%m%y%H%M%S"),
        STAMPY_VERSION="stampy-1.0.27",
        PICARD_VERSION="picard-tools-1.128",
        FAST="--end-to-end --sensitive",
        SLOW="--local --sensitive-local",
        BAM_PREFIX=align_data.bam_prefix(),

        TRIMMED_FQ1=trimmed_fq[0],
        TRIMMED_FQ2=trimmed_fq[1],
        DUPEMARK_BAM=align_data.dupl_mark_bam(),
        DUPEMARK_FLAGSTAT=align_data.dupl_mark_bam_flagstat(),
        DUPEMARK_BED=align_data.dupl_mark_bed(),

        STAMPY_TMP_BAM=align_data.tmp_bam(),
        STAMPY_DUPEMARK_CLEANSAM_TMP_BAM=align_data.dupl_mark_tmp_cleansam_bam(),
        STAMPY_DUPEMARK_TMP_BAM=align_data.dupl_mark_tmp_bam(),

        DUPEMARK_REALN_BAM=realn_data.dupl_mark_realn_bam(),
        DUPEMARK_REALN_FLAGSTAT=realn_data.dupl_mark_realn_bam_flagstat(),
        DUPEMARK_REALN_BED=realn_data.dupl_mark_realn_bed(),
        DUPEMARK_BAM_FOR_INDER_REALN_INTERVALS=realn_data.reports_path(
            realn_data.bam_prefix() + ".dupemk.bam.ForIndelRealigner.intervals"),

        DISCORDANT_SAM=align_data.discordant_sam(),
        DISCORDANT_BAM=align_data.discordant_bam(),
        SPLITREAD_SAM=align_data.splitread_sam(),
        SPLITREAD_BAM=align_data.splitread_bam(),
        UNMAPPED_FASTQ=align_data.unmapped_fastq(),

        TMP_DIR=align_data.tmp_dir(),

        REFFASTA=genomebuild.ref_fasta(),
        KNOWN_INDELS=genomebuild.known_indels(),
        REFFASTA=genomebuild.ref_fasta(),
        DBSNP_RECAB=genomebuild.dbsnp_recab(),
        KNOWN_SNPS_b138=genomebuild.known_snps_b138(),
        KNOWN_SNPS_OMNI=genomebuild.known_snps_omni(),
        KNOWN_SNPS_1000G=genomebuild.known_snps_1000g(),

        BSQR_BAM_IN=bsqr_data.bsqr_bam_in(),
        BSQR_BAM_OUT=bsqr_data.bsqr_bam_out(),
        RECAL_DATA_TABLE=bsqr_data.reports_path(bsqr_data.bam_prefix() + ".recal_data.table"),

        VC_BAM_IN=vc_data.vc_bam_in(),
        VC_FILTERED_BAM=vc_data.vc_filtered_bam(),
        MAPPED_READS_BED=vc_data.reports_path(os.path.basename(vc_data.vc_bam_in()) + ".mapped.reads.bed"),
        GENOMECOV_BED=vc_data.reports_path(os.path.basename(vc_data.vc_bam_in()) + ".genomecov.bed"),
        BASE_QUAL="20",
        MAP_QUAL="20",
        COVERAGE_MIN="2",
        RAW_VCF=vc_data.raw_vcf(),
        RAW_VCF_GZ=vc_data.raw_vcf_gz(),
        VCF_GZ=vc_data.vcf_gz(),
        FREEBAYES_RAW_VCF=vc_data.raw_vcf("freebayes-parallel"),
        FREEBAYES_RAW_VCF_GZ=vc_data.raw_vcf_gz("freebayes-parallel"),
        FREEBAYES_VCF_GZ=vc_data.vcf_gz("freebayes-parallel"),
        HAPLOTYPE_CALLER_RAW_VCF=vc_data.raw_vcf("HaplotypeCaller"),
        HAPLOTYPE_CALLER_VCF_GZ=vc_data.raw_vcf_gz("HaplotypeCaller"),
        HAPLOTYPE_CALLER_VCF_GZ=vc_data.vcf_gz("HaplotypeCaller")
    )

    dict.set_lazy("PLATFORM_UNIT", find_platform_unit, fqc_data.fastq_files()[0])
    return dict


def select_genome(row, projects_home):
    genome = genome_build.select(row.genomebuild(), projects_home)
    if genome is None:
        raise ValueError("No genome selected. Choose one of [b37] or [hg19]")
    return genome


def find_platform_unit(fastq_file):
    cmd = "zcat %s | head -1 | sed 's/:/\\t/' - | cut -f 1 | sed 's/@//g' - " % fastq_file
    logger().debug("platform_info=[%s]" % cmd)
    out = subprocess.check_output(cmd, shell=True, preexec_fn=lambda: signal(SIGPIPE, SIG_DFL))
    return out.strip()
