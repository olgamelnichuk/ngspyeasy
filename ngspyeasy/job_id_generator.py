#!/usr/bin/env python

import itertools
import datetime

COUNTER = itertools.count()

TIMESTAMP = datetime.datetime.now().strftime("%d%m%y%H%M%S")


def get_next(tokens):
    return  "_".join([TIMESTAMP, str(COUNTER.next())] + tokens)