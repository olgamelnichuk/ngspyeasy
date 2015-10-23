#!/usr/bin/env python
import sys
import logging

import os.path

FORMATTER = logging.Formatter('%(asctime)s %(threadName)s %(module)s %(levelname)s: %(message)s')


def init_logger(logfile, verbose):
    log_dir = os.path.dirname(logfile)
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(file_handler(logfile))
    logger.addHandler(console_handler())
    return logger


def file_handler(logfile):
    handler = logging.FileHandler(logfile)
    handler.setFormatter(FORMATTER)
    return handler


def console_handler():
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(FORMATTER)
    return handler


def logger():
    logger = logging.getLogger()
    if len(logger.handlers) == 0:
        logger.addHandler(console_handler())
        logger.setLevel(logging.DEBUG)
    return logger