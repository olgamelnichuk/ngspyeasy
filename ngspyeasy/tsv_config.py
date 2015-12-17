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

import csv

import os.path


def parse(tsv_path):
    if not os.path.exists(tsv_path):
        raise IOError("File %s does not exist" % tsv_path)

    filename = os.path.basename(tsv_path)

    if os.stat(tsv_path).st_size == 0:
        return TsvConfig([])

    rows = []
    with open(tsv_path, 'r') as tsv:
        reader = csv.reader(utf_8_encoder(tsv), delimiter="\t")
        try:
            for row in reader:
                rows.append(row)
        except csv.Error as e:
            raise ValueError('TSV Format Error: file %s, line %d: %s' % (tsv_path, reader.line_num, e))

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

        header = rows[0]
        self.rows = [TsvConfigRow("sample_" + str(i), x, header) for i, x in enumerate(rows[1:])]

    def all_rows(self):
        for r in self.rows:
            yield r

    def is_empty(self):
        return len(self.rows) == 0

    def row_size(self):
        return len(self.rows)

    def col_size(self):
        return 0 if self.is_empty() else self.row_at(0).size()

    def row_at(self, index):
        return self.rows[index] if 0 <= index <= self.row_size() else None


class TsvConfigRow:
    def __init__(self, id, values, headers):
        d = dict()
        for i, h in enumerate(headers, start=0):
            d[h.lower()] = values[i]

        self.__dict__.update(d)
        self.__id__ = id

    def size(self):
        return len(self.__dict__)

    def __str__(self):
        return str(self.__dict__)

    def id(self):
        return self.__id__
