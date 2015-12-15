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

import sys
import logging
import datetime

import os.path

FORMATTER = logging.Formatter('%(asctime)s %(threadName)s %(module)s %(levelname)s: %(message)s')


def init_main_logger(log_dir):
    logname = "%s_ngseasy.log" % datetime.datetime.now().strftime("%d%m%y")
    init_logger(os.path.join(log_dir, logname), verbose=True)


def init_sample_logger(log_dir, sample_id):
    logname = "%s_ngseasy@%s.log" % (datetime.datetime.now().strftime("%d%m%y"), sample_id)
    log_path = os.path.join(log_dir, logname)
    init_logger(log_path, verbose=True)


def init_logger(logfile, verbose=False):
    log_dir = os.path.dirname(logfile)
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger("main")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(file_handler(logfile))
    logger.addHandler(console_handler())

    file_only = logging.getLogger("file-only")
    file_only.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_only.addHandler(file_handler(logfile))

    return logger


def file_handler(logfile):
    handler = logging.FileHandler(logfile)
    handler.setFormatter(FORMATTER)
    return handler


def console_handler():
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(FORMATTER)
    return handler


def logger(file_only=False):
    if file_only:
        return logging.getLogger("file-only")

    logger = logging.getLogger("main")
    if len(logger.handlers) == 0:
        logger.addHandler(console_handler())
        logger.setLevel(logging.DEBUG)
    return logger
