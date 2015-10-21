#!/usr/bin/env python

import re
import projects_dir
from utils import Bunch

illumina_fastq = [r'(.*_L.*_R[1,2]_[0-9][0-9][0-9])\.fastq\.gz',
                  r'(.*_L.*_R[1,2]_[0-9][0-9][0-9][0-9])\.fastq\.gz']

other_fastq = [r'(.*)_([1,2])\.fastq\.gz',
               r'(.*)_R([1,2])\.fastq\.gz',
               r'(.*)_([1,2])\.fq\.gz',
               r'(.*)_R([1,2])\.fq\.gz']


def matched_fastq_name(filename):
    for p in illumina_fastq:
        match = re.match(p, filename)
        if match:
            return match.group(1), None, "illumina"

    for p in other_fastq:
        match = re.match(p, filename)
        if match:
            return match.group(1), match.group(2), "other"

    return None, None, None


def parsed_fastq_name(fastq_name):
    (prefix, num, fq_type) = matched_fastq_name(fastq_name)
    if prefix is None:
        return None

    name = [x for x in [prefix, num] if x is not None]
    return Bunch(name=name, type=fq_type)


def create(row, projects_home):
    sample_dir = projects_home.sample_dir(row.project_id(), row.sample_id())

    for fastq_name in [row.fastq1(), row.fastq2()]:
        if parsed_fastq_name(fastq_name) is None:
            raise ValueError("FastQ naming format is not recognised for name: %s", fastq_name)

    return row, sample_dir


def fastqc_data(row, projects_home):
    return FastQCData(*create(row, projects_home))


def trimmomatic_data(row, projects_home):
    return TrimmomaticData(*create(row, projects_home))


def alignment_data(row, projects_home):
    return AlignmentData(*create(row, projects_home))


def realn_data(row, projects_home):
    return RealnData(*create(row, projects_home))


def bsqr_data(row, projects_home):
    return BsqrData(*create(row, projects_home))


def vc_data(row, projects_home):
    return VarCallerData(*create(row, projects_home))


class SampleData(projects_dir.SampleDir):
    def __init__(self, row, sample_dir):
        super(SampleData, self).__init__(sample_dir)
        self._row = row

    def fastq_names(self):
        return [self._row.fastq1(), self._row.fastq2()]

    def row(self):
        return self._row


class FastQCData(SampleData):
    def fastq_files(self):
        return [self.fastq_path(x) for x in self.fastq_names()]

    def fastqc_htmls(self):
        parsed_names = [parsed_fastq_name(x) for x in self.fastq_names()]
        htmls = ["_".join(x.name + ["fastqc.html"]) for x in parsed_names]
        return [self.fastq_path(x) for x in htmls]


class TrimmomaticData(FastQCData):
    def paired_fastq(self):
        return self.result_fastq("filtered.fastq.gz")

    def unpaired_fastq(self):
        return self.result_fastq("unpaired.fastq.gz")

    def result_fastq(self, suffix):
        r = self.row()
        results = [".".join(
            [r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.trim() + "_" + str(x), suffix])
                   for x in [1, 2]]
        return [self.fastq_path(x) for x in results]

    def trim_fastqc_htmls(self):
        return [x.replace("fastq.gz", "fastqc.html") for x in self.unpaired_fastq() + self.paired_fastq()]


class AlignmentData(TrimmomaticData):
    def notrim(self):
        return self.row().trim() == "no-trim"

    def trimmed_fastq_files(self):
        return self.fastq_files() if self.notrim() else self.paired_fastq()

    def bam_prefix(self):
        r = self.row()
        return ".".join([r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.ngs_platform(), r.trim(), r.aligner(),
                         r.genomebuild()])

    def discordant_sam(self):
        return self.alignments_path(self.bam_prefix() + ".discordant.sam")

    def discordant_bam(self):
        return self.alignments_path(self.bam_prefix() + ".discordant.bam")

    def splitread_sam(self):
        return self.alignments_path(self.bam_prefix() + ".splitread.sam")

    def splitread_bam(self):
        return self.alignments_path(self.bam_prefix() + ".splitread.bam")

    def unmapped_fastq(self):
        return self.alignments_path(self.bam_prefix() + ".unmapped.fastq")

    def dupl_mark_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.bam")

    def dupl_mark_bam_flagstat(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.bam.flagstat")

    def dupl_mark_bed(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.bed")

    def tmp_bam(self):
        return self.alignments_path(self.bam_prefix() + ".tmp.bam")

    def dupl_mark_tmp_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.tmp.bam")

    def dupl_mark_tmp_cleansam_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.tmpcleansam.bam")


class RealnData(AlignmentData):
    def dupl_mark_realn_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.bam" % self.row().realn())

    def dupl_mark_realn_bam_flagstat(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.%s.bam.flagstat" % self.row().realn())

    def dupl_mark_realn_bed(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.%s.bed" % self.row().realn())


class BsqrData(RealnData):
    def bsqr_bam_in(self):
        r = self.row()
        if (r.bsqr() == "bam-bsqr" and r.realn() == "bam-realn") or (
                        r.bsqr() == "gatk-bsqr" and r.realn() == "gatk-realn"):
            return self.dupl_mark_realn_bam()
        return self.dupl_mark_bam()

    def bsqr_bam_out(self):
        r = self.row()
        if (r.bsqr() == "bam-bsqr" and r.realn() == "bam-realn") or (
                        r.bsqr() == "gatk-bsqr" and r.realn() == "gatk-realn"):
            return self.dupl_mark_realn_bsqr_bam()
        return self.dupl_mark_bsqr_bam()

    def dupl_mark_realn_bsqr_bam(self):
        r = self.row()
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.%s.bam" % (r.realn(), r.bsqr()))

    def dupl_mark_bsqr_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.bam" % self.row().bsqr())


class VarCallerData(BsqrData):
    def vc_bam_in(self):
        realn = self.row().realn()
        bsqr = self.row().bsqr()
        if realn == "no-realn" and bsqr == "no-bsqr":
            return self.dupl_mark_bam()
        elif realn == "no-bsqr":
            return self.dupl_mark_realn_bam()
        else:
            return self.bsqr_bam_out()

    def vc_filtered_bam(self):
        return self.bsqr_bam_in().replace(".bam", ".filtered.bam")

    def raw_vcf(self, varcaller=None):
        if varcaller is None:
            varcaller = self.row().varcaller()
        return self.vcf_path(self.bam_prefix() + ".raw.snps.indels.%s.vcf" % varcaller)

    def raw_vcf_gz(self, varcaller=None):
        return self.raw_vcf(varcaller) + ".gz"

    def vcf_gz(self, varcaller=None):
        if varcaller is None:
            varcaller = self.row().varcaller()

        return self.vcf_path(self.bam_prefix() + ".snps.indels.%s.vcf.gz" % varcaller)
