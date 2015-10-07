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


def parse_fastq_name(filename):
    for p in illumina_fastq:
        match = re.match(p, filename)
        if match:
            return match.group(1), None, "illumina"

    for p in other_fastq:
        match = re.match(p, filename)
        if match:
            return match.group(1), match.group(2), "other"

    return None, None, None


def normalize_fastq(path):
    (dir, name) = os.path.split(path)
    (prefix, num, fq_type) = parse_fastq_name(name)

    if prefix is None:
        raise ValueError("Fastq naming format not recognised for file: %s", path)

    name = [x for x in [prefix, num] if x is not None]

    return Bunch(path=path, type=fq_type, name=name)


def create(row, projects_home):
    sample_dir = projects_dir.sample_dir(projects_home, row.project_id(), row.sample_id())

    fastq = [row.fastq1(), row.fastq2()]
    fastq = [projects_dir.fastq_path(sample_dir, x) for x in fastq]

    for fastq_file in fastq:
        if not os.path.isfile(fastq_file):
            raise IOError("FastQ file not found: %s", fastq_file)

    fastq_normalized = map(lambda x: normalize_fastq(x), fastq)
    fastq_types = set(map(lambda x: x.type, fastq_normalized))

    if len(fastq_types) > 1:
        raise ValueError("FastQ file formats are not the same: %s" % str(fastq_types))

    return SampleData(row, sample_dir, fastq_normalized)


class SampleData(object):
    def __init__(self, row, sample_dir, fastq_normalized):
        self.sample_dir = sample_dir
        self.fastq = fastq_normalized
        self.row = row

    def tmp_dir(self):
        return projects_dir.sample_tmp_dir(self.sample_dir)

    def fastq_dir(self):
        return projects_dir.sample_fastq_dir(self.sample_dir)

    def fastq_files(self):
        return [x.path for x in self.fastq]

    def fastqc_results(self):
        results = ["_".join(x.name + ["fastqc.html"]) for x in self.fastq]
        return [os.path.join(self.fastq_dir(), x) for x in results]

    def trimmomatic_paired_results(self):
        return self.trimmomatic_results("filtered.fastq.gz")

    def trimmomatic_unpaired_results(self):
        return self.trimmomatic_results("unpaired.fastq.gz")

    def trimmomatic_results(self, suffix):
        r = self.row
        results = [".".join(
            [r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.trim() + "_" + str(x), suffix])
                   for x in [1, 2]]
        return [os.path.join(self.fastq_dir(), x) for x in results]

    def bam_prefix(self):
        r = self.row
        return ".".join([r.sample_id(), r.ngs_type(), r.dna_prep_library_id(), r.ngs_platform(), r.trim(), r.aligner(),
                         r.genomebuild()])

    def alignments_path(self, filename):
        dir = projects_dir.sample_alignments_dir(self.sample_dir)
        return os.path.join(dir, filename)
