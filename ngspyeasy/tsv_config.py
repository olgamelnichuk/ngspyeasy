#!/usr/bin/env python

import csv

def parse(tsvfile, log):
    rows = []
    with open(tsvfile) as tsv:
        reader = csv.reader(tsv, delimeter="\t")
        try:
            for row in reader:
                rows.append(row)
        except csv.Error as e:
            log.error('file %s, line %d: %s' % (tsvfile, reader.line_num, e))
            return None
    return TsvConfig(rows)


class TsvConfig:

    def __init__(self, rows):
        self.rows = rows

    def isEmpty(self):
        len(self.rows) == 0