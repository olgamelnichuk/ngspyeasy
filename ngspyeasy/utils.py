#!/usr/bin/env python
import os

import re

illumina_fastq = [r'(.*_L.*_R[1,2]_[0-9][0-9][0-9])\.fastq\.gz',
                  r'(.*_L.*_R[1,2]_[0-9][0-9][0-9][0-9])\.fastq\.gz']

other_fastq = [r'(.*)_([1,2])\.fastq\.gz',
               r'(.*)_R([1,2])\.fastq\.gz',
               r'(.*)_([1,2])\.fq\.gz',
               r'(.*)_R([1,2])\.fq\.gz']


def uniq_set(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def recognize_fastq_naming(filename):
    for p in illumina_fastq:
        match = re.match(p, filename)
        if match:
            return match.group(1), None, "illumina"

    for p in other_fastq:
        match = re.match(p, filename)
        if match:
            return match.group(1), match.group(2), "other"

    return None, None, None


def recognize_fastq(path):
    (dir, name) = os.path.split(path)
    (prefix, num, fq_type) = recognize_fastq_naming(name)

    if prefix is None:
        raise ValueError("Fastq naming format not recognised for file: %s", path)

    result = [prefix, num, "fastqc.html"]
    result = [x for x in result if x is not None]

    return Bunch(type=fq_type, result=os.path.join(dir, "_".join(result)))


class Enum(object):
    def __init__(self, *keys):
        self.__dict__.update(zip(keys, range(len(keys))))


class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__