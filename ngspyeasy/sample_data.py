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

    def fastqc_output(self):
        parsed_names = [parsed_fastq_name(x) for x in self.fastq_names()]
        htmls = ["_".join(x.name + ["fastqc.html"]) for x in parsed_names]
        return [self.dir.fastq_path(x) for x in htmls]

    def trimmomatic_paired_output(self):
        return self.trimmomatic_output("filtered.fastq.gz")

    def trimmomatic_unpaired_output(self):
        return self.trimmomatic_output("unpaired.fastq.gz")

    def trimmomatic_output(self, suffix):
        r = self.row
        results = [".".join(
            [r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.trim() + "_" + str(x), suffix])
                   for x in [1, 2]]
        return [self.fastq_path(x) for x in results]

    def bam_prefix(self):
        r = self.row
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

    def dupl_mark_filtered_bam(self):
        return self.alignments_path(self.bam_prefix() + ".dupemk.filtered.bam")

    def dupl_mark_bam_flagstat(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.bam.flagstat")

    def dupl_mark_bed(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.bed")

    def dupl_mark_realn_bam(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.bam" % realn)

    def dupl_mark_realn_bsqr_bam(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.%s.bam" % (realn, self.row.bsqr()))

    def dupl_mark_realn_bsqr_filtered_bam(self, realn=None):
        realn = self.row.realn() if realn is None else realn
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.%s.filtered.bam" % (realn, self.row.bsqr()))

    def dupl_mark_bsqr_bam(self, bsqr=None):
        bsqr = self.row.bsqr() if bsqr is None else bsqr
        return self.alignments_path(self.bam_prefix() + ".dupemk.%s.bam" % bsqr)

    def dupl_mark_realn_bam_flagstat(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.%s.bam.flagstat" % self.row.realn())

    def dupl_mark_realn_bed(self):
        return self.reports_path(self.bam_prefix() + ".dupemk.%s.bed" % self.row.realn())
