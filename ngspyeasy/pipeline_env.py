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
import inspect
from utils import LazyDict


def as_dict(row, projects_home):
    return environment_dict(PipelineEnvironment(row, projects_home))


def as_test_dict():
    return environment_dict(None)


def environment_dict(obj):
    dict = LazyDict()
    for (name, method) in [(n, m) for (n, m) in inspect.getmembers(PipelineEnvironment) if non_inherited_method(n, m)]:
        dict.set_lazy(name.upper(), "DUMMY" if obj is None else getattr(obj, name))
    return dict


def non_inherited_method(method_name, method_obj):
    return inspect.ismethod(method_obj) and not hasattr(PipelineEnvironmentBase, method_name)


class PipelineEnvironmentBase(object):
    def select_genome(self, row, projects_home):
        genome = genome_build.select(row.genomebuild(), projects_home)
        if genome is None:
            raise ValueError("No genome selected. Choose one of [b37] or [hg19]")
        return genome

    def find_platform_unit(self, fastq_file):
        cmd = "zcat %s | head -1 | sed 's/:/\\t/' - | cut -f 1 | sed 's/@//g' - " % fastq_file
        logger().debug("platform_info=[%s]" % cmd)
        out = subprocess.check_output(cmd, shell=True, preexec_fn=lambda: signal(SIGPIPE, SIG_DFL))
        return out.strip()


class PipelineEnvironment(PipelineEnvironmentBase):
    def __init__(self, row, projects_home):
        self.row = row
        self.fqc_data = sample.fastqc_data(row, projects_home)
        self.trim_data = sample.trimmomatic_data(row, projects_home)
        self.align_data = sample.alignment_data(row, projects_home)
        self.realn_data = sample.realn_data(row, projects_home)
        self.bsqr_data = sample.bsqr_data(row, projects_home)
        self.vc_data = sample.vc_data(row, projects_home)

        self.fq = self.fqc_data.fastq_files()
        self.trimmed_fq = self.align_data.trimmed_fastq_files()
        self.genomebuild = self.select_genome(row, projects_home)

    def ncpu(self):
        return self.row.ncpu()

    def ngs_platform(self):
        return self.row.ngs_platform()

    def dna_prep_library_id(self):
        return self.row.dna_prep_library_id()

    def rundate(self):
        return datetime.datetime.now().strftime("%d%m%y%H%M%S")

    def stampy_version(self):
        return "stampy-1.0.27"

    def picard_version(self):
        return "picard-tools-1.128"

    def fast(self):
        return "--end-to-end --sensitive"

    def slow(self):
        return "--local --sensitive-local"

    def bam_prefix(self):
        return self.align_data.bam_prefix()

    def trimmed_fq1(self):
        return self.trimmed_fq[0]

    def trimmed_fq2(self):
        return self.trimmed_fq[1]

    def dupemark_bam(self):
        return self.align_data.dupl_mark_bam()

    def dupemark_flagstat(self):
        return self.align_data.dupl_mark_bam_flagstat()

    def dupemark_bed(self):
        return self.align_data.dupl_mark_bed()

    def stampy_tmp_bam(self):
        return self.align_data.tmp_bam()

    def stampy_dupemark_tmp_bam(self):
        return self.align_data.dupl_mark_tmp_bam()

    def stampy_dupemark_cleansam_tmp_bam(self):
        return self.align_data.dupl_mark_tmp_cleansam_bam()

    def dupemark_realn_bam(self):
        return self.realn_data.dupl_mark_realn_bam()

    def dupemark_realn_flagstat(self):
        return self.realn_data.dupl_mark_realn_bam_flagstat()

    def dupemark_realn_bed(self):
        return self.realn_data.dupl_mark_realn_bed()

    def dupemark_bam_for_inder_realn_intervals(self):
        return self.realn_data.reports_path(
            self.realn_data.bam_prefix() + ".dupemk.bam.ForIndelRealigner.intervals")

    def discordant_sam(self):
        return self.align_data.discordant_sam()

    def discordant_bam(self):
        return self.align_data.discordant_bam()

    def splitread_sam(self):
        return self.align_data.splitread_sam()

    def splitread_bam(self):
        return self.align_data.splitread_bam()

    def unmapped_fastq(self):
        return self.align_data.unmapped_fastq()

    def tmp_dir(self):
        return self.align_data.tmp_dir()

    def reffasta(self):
        return self.genomebuild.ref_fasta()

    def known_indels(self):
        return self.genomebuild.known_indels()

    def dbsnp_recab(self):
        return self.genomebuild.dbsnp_recab()

    def known_snps_b138(self):
        return self.genomebuild.known_snps_b138()

    def known_snps_omni(self):
        return self.genomebuild.known_snps_omni()

    def known_snps_1000G(self):
        return self.genomebuild.known_snps_1000g()

    def genomeindex(self):
        return self.genomebuild.genome_index()

    def novoindex(self):
        return self.genomebuild.novoindex()

    def bsqr_bam_in(self):
        return self.bsqr_data.bsqr_bam_in()

    def bsqr_bam_out(self):
        return self.bsqr_data.bsqr_bam_out()

    def recal_data_table(self):
        return self.bsqr_data.reports_path(self.bsqr_data.bam_prefix() + ".recal_data.table"),

    def vc_bam_in(self):
        return self.vc_data.vc_bam_in()

    def vc_filtered_bam(self):
        return self.vc_data.vc_filtered_bam()

    def mapped_reads_bed(self):
        return self.vc_data.reports_path(os.path.basename(self.vc_data.vc_bam_in()) + ".mapped.reads.bed")

    def genomecov_bed(self):
        return self.vc_data.reports_path(os.path.basename(self.vc_data.vc_bam_in()) + ".genomecov.bed")

    def base_qual(self):
        return "20"

    def map_qual(self):
        return "20"

    def coverage_min(self):
        return "2"

    def raw_vcf(self):
        return self.vc_data.raw_vcf()

    def raw_vcf_gz(self):
        return self.vc_data.raw_vcf_gz()

    def vcf_gz(self):
        return self.vc_data.vcf_gz()

    def freebayes_raw_vcf(self):
        return self.vc_data.raw_vcf("freebayes-parallel")

    def freebayes_raw_vcf_gz(self):
        return self.vc_data.raw_vcf_gz("freebayes-parallel")

    def freebayes_vcf_gz(self):
        return self.vc_data.vcf_gz("freebayes-parallel")

    def haplotype_caller_raw_vcf(self):
        return self.vc_data.raw_vcf("HaplotypeCaller")

    def haplotype_caller_raw_vcf_gz(self):
        return self.vc_data.raw_vcf_gz("HaplotypeCaller")

    def haplotype_caller_vcf_gz(self):
        return self.vc_data.vcf_gz("HaplotypeCaller")

    def platypus_raw_vcf(self):
        return self.vc_data.raw_vcf("platypus")

    def platypus_raw_vcf_gz(self):
        return self.vc_data.raw_vcf_gz("platypus")

    def platypus_vcf_gz(self):
        return self.vc_data.vcf_gz("platypus")

    def platform_unit(self):
        return self.find_platform_unit(self.fqc_data.fastq_files()[0])

    def k_stats(self):
        return self.align_data.alignments_path(self.align_data.bam_prefix() + ".K.stats")

    def refdir(self):
        return self.genomebuild.refdir()

    def fastq_dir(self):
        return self.fqc_data.fastq_dir()

    def fq1(self):
        return self.fqc_data.fastq_files()[0]

    def fq2(self):
        return self.fqc_data.fastq_files()[1]

    def paired_fq1(self):
        return self.trim_data.paired_fastq()[0]

    def paired_fq2(self):
        return self.trim_data.paired_fastq()[1]

    def unpaired_fq1(self):
        return self.trim_data.unpaired_fastq()[0]

    def unpaired_fq2(self):
        return self.trim_data.unpaired_fastq()[1]

    def adapter_fa(self):
        return self.genomebuild.adapter_fa()

    def gtmodegatk(self):
        return "EMIT_VARIANTS_ONLY"

    def filter_duplicates(self):
        ngs_type = self.row.ngs_type()
        return 1 if ngs_type == "TGS" or ngs_type == "WEX" else 0