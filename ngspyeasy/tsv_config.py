#!/usr/bin/env python

import csv
import os.path
from utils import Enum
from logger import log_warn, log_error

columns = Enum('PROJECT_ID', 'SAMPLE_ID', 'FASTQ1', 'FASTQ2', 'PROJECT_DIR', 'DNA_PREP_LIBRARY_ID', 'NGS_PLATFORM',
               'NGS_TYPE', 'BAIT', 'CAPTURE', 'GENOMEBUILD', 'FASTQC', 'TRIM', 'REALN', 'BSQR', 'ALIGNER', 'VARCALLER',
               'CNV', 'ANNOTATOR', 'CLEANUP', 'NCPU', 'VERSION', 'NGSUSER')


def parse(tsvfile):
    if not os.path.exists(tsvfile):
        log_error("[TSV parser]: file %s doesn't exist", tsvfile)
        return None

    if os.stat(tsvfile).st_size == 0:
        log_warn("[TSV parser]: file %s is empty", tsvfile)
        return TsvConfig([])

    rows = []
    with open(tsvfile, 'r') as tsv:
        reader = csv.reader(utf_8_encoder(tsv), delimiter="\t")
        try:
            for row in reader:
                rows.append(row)
        except csv.Error as e:
            log_error('file %s, line %d: %s' % (tsvfile, reader.line_num, e))
            return None
    return TsvConfig(rows)


def utf_8_encoder(unicode_csv_data):
    for line in unicode_csv_data:
        yield line.encode('utf-8')


class TsvConfig:
    def __init__(self, rows):
        self.rows = []
        self.header = None

        if (rows is None) or len(rows) == 0:
            return

        header = TsvConfigRow(rows[0])
        data_rows = rows[0:]

        if header.get_project_id() == 'PROJECT_ID':
            self.header = header
            data_rows = rows[1:]

        self.rows = map(lambda x: TsvConfigRow(x), data_rows)

    def has_header(self):
        return self.header is not None

    def all_rows(self):
        return self.rows[0:]

    def is_empty(self):
        return len(self.rows) == 0

    def row_size(self):
        return len(self.rows)

    def col_size(self):
        return 0 if self.is_empty() else self.row_at(0).size()

    def row_at(self, index):
        return self.rows[index] if 0 <= index <= self.row_size() else None


class TsvConfigRow:
    def __init__(self, values):
        self.values = values

    def size(self):
        return len(self.values)

    def value_at(self, index):
        return self.values[index] if 0 <= index <= len(self.values) else None

    def get_project_id(self):
        return self.value_at(columns.PROJECT_ID)

    def get_sample_id(self):
        return self.value_at(columns.SAMPLE_ID)
