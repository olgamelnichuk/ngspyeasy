#!/usr/bin/env python

import os
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
    sample_dir = projects_dir.SampleDir(
        projects_home.sample_dir(row.project_id(), row.sample_id()))

    for fastq_name in [row.fastq1(), row.fastq2()]:
        if parsed_fastq_name(fastq_name) is None:
            raise ValueError("FastQ naming format is not recognised for name: %s", fastq_name)

    sample = SampleData(row, sample_dir)

    for fastq_file in sample.fastq_files():
        if not os.path.isfile(fastq_file):
            raise IOError("FastQ file not found: %s", fastq_file)

    return sample


class SampleData(object):
    def __init__(self, row, sample_dir):
        self.dir = sample_dir
        self.row = row

    def tmp_dir(self):
        return self.dir.tmp_dir()

    def fastq_dir(self):
        return self.dir.fastq_dir()

    def alignments_dir(self):
        return self.dir.alignments_dir()

    def fastq_path(self, filename):
        return os.path.join(self.dir.fastq_dir(), filename)

    def alignments_path(self, filename):
        return os.path.join(self.dir.alignments_dir(), filename)

    def reports_path(self, filename):
        return os.path.join(self.dir.reports_dir(), filename)

    def fastq_names(self):
        return [self.row.fastq1(), self.row.fastq2()]

    def fastq_files(self):
        return [self.dir.fastq_path(x) for x in self.fastq_names()]

    def bam_prefix(self):
        r = self.row
        return ".".join([r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.ngs_platform(), r.trim(), r.aligner(),
                         r.genomebuild()])

    def fastq_data(self):
        return FastQCData(self.row, self.dir)

    def trimmomatic_data(self):
        return TrimmomaticData(self.row, self.dir)

    def alignment_data(self):
        return AlignmentData(self.row, self.dir)

    def realn_data(self):
        return RealnData(self.row, self.dir)

    def bsqr_data(self):
        return BsqrData(self.row, self.dir)

    def vc_data(self):
        return VarCallerData(self.row, self.dir)


class FastQCData(SampleData):
    def fastqc_htmls(self):
        parsed_names = [parsed_fastq_name(x) for x in self.fastq_names()]
        htmls = ["_".join(x.name + ["fastqc.html"]) for x in parsed_names]
        return [self.fastq_path(x) for x in htmls]


class TrimmomaticData(SampleData):
    def paired_fastq(self):
        return self.result_fastq("filtered.fastq.gz")

    def unpaired_fastq(self):
        return self.result_fastq("unpaired.fastq.gz")

    def result_fastq(self, suffix):
        r = self.row
        results = [".".join(
            [r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.trim() + "_" + str(x), suffix])
                   for x in [1, 2]]
        return [self.fastq_path(x) for x in results]


class AlignmentData(SampleData):
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


class RealnData(SampleData):
    def dupl_mark_realn_bam(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.bam" % realn)

    def dupl_mark_realn_bam_flagstat(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.reports_path(self.bam_prefix() + ".dupemk.%s.bam.flagstat" % realn)

    def dupl_mark_realn_bed(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.reports_path(self.bam_prefix() + ".dupemk.%s.bed" % realn)


class BsqrData(SampleData):
    def dupl_mark_realn_bsqr_bam(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.%s.bam" % (realn, self.row.bsqr()))

    def dupl_mark_bsqr_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.bam" % self.row.bsqr())


class VarCallerData(SampleData):
    def dupl_mark_filtered_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.filtered.bam")

    def dupl_mark_realn_bsqr_filtered_bam(self):
        return self.alignments_path(
            self.bam_prefix() + ".dupemk.%s.%s.filtered.bam" % (self.row.realn(), self.row.bsqr()))

    def bam_file(self):
        if self.row.realn() == "no-realn" and self.row().bsqr() == "no-bsqr":
            return self.alignment_data().dupl_mark_bam()
        return self.bsqr_data().dupl_mark_realn_bsqr_bam()

    def filtered_bam_file(self):
        if self.row.realn() == "no-realn" and self.row().bsqr() == "no-bsqr":
            return self.dupl_mark_filtered_bam()
        return self.dupl_mark_realn_bsqr_filtered_bam()
