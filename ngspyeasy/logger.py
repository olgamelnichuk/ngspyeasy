#!/usr/bin/env python
import sys
import logging
import os.path

LOGGER_NAME = "ngspyeasy"

def init_logger(logfile, verbose):
    log_dir = os.path.dirname(logfile)
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger(LOGGER_NAME)
    formatter = logging.Formatter('%(asctime)s %(threadName)s %(name)s %(levelname)s: %(message)s')

    file_handler = logging.FileHandler(logfile)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def get_logger():
    return logging.getLogger(LOGGER_NAME)
