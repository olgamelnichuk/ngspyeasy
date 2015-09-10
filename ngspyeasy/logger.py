#!/usr/bin/env python

import os.path
import logging

def create(logdir, tsv_name, verbose):
    version = "0.1"
    rundate = "RUNDATE"
    logfile = os.path.join(logdir, "ngspyeasy." + version + "." + tsv_name + "." + rundate)
    return Logger(logfile, verbose)


class Logger:

    def __init__(self, filename=None, verbose=False):
        logger = logging.getLogger()
        formatter = logging.Formatter('%(asctime)s : %(message)s')
        fileHandler = logging.FileHandler(filename, mode='w')
        fileHandler.setFormatter(formatter)
        streamHandler = logging.StreamHandler()
        streamHandler.setFormatter(formatter)

        logger.setLevel()
        logger.addHandler(fileHandler)
        logger.addHandler(streamHandler)
        self.logger = logger

    def error(self, msg):
        self.logger.error(msg)

    def debug(self, msg):
        self.logger.debug(msg)

    def info(self, msg):
        self.logger.info(msg)