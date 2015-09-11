#!/usr/bin/env python

import csv
import os.path
from logger import log_warn, log_error


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
        self.rows = map(lambda x: TsvConfigRow(x), rows)

    def is_empty(self):
        return len(self.rows) == 0

    def row_size(self):
        return len(self.rows)

    def col_size(self):
        return 0 if self.is_empty() else self.get(0).size()

    def get(self, index):
        return self.rows[index] if 0 <= index <= self.row_size() else None


class TsvConfigRow:
    def __init__(self, values):
        self.values = values

    def size(self):
        return len(self.values)
