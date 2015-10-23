#!/usr/bin/env python
import sys
import logging

import os.path

FORMATTER = logging.Formatter('%(asctime)s %(threadName)s %(module)s %(levelname)s: %(message)s')


def logger_name(name=None):
    return "ngspyeasy" + ("" if name is None else "." + name)


def init_logger(logfile, verbose, name=None):
    log_dir = os.path.dirname(logfile)
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger(logger_name(name))
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


def get_logger(name=None):
    logger = logging.getLogger(logger_name(name))
    if len(logger.handlers) == 0:
        logger.addHandler(console_handler())
        logger.setLevel(logging.DEBUG)
    return logger


def log_debug(msg):
    get_logger().debug(msg)


def log_info(msg):
    get_logger().info(msg)


def log_error(msg):
    get_logger().error(msg)


def log_exception(ex):
    get_logger().exception(ex)
