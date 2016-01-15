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
import multiprocessing
import threading
import traceback

import os.path

FORMATTER = logging.Formatter('%(asctime)s %(threadName)s %(module)s %(levelname)s: %(message)s')


def init_main_logger(log_dir, verbose=True):
    logname = "%s_ngseasy.log" % datetime.datetime.now().strftime("%d%m%y")
    logfile = os.path.join(log_dir, logname)

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger("root")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(multi_proc_handler(logfile))
    logger.addHandler(console_handler())

    file_only = logging.getLogger("file-only")
    file_only.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_only.addHandler(multi_proc_handler(logfile))
    return logger


def init_play_run_logger(log_dir, run_id, verbose=True):
    logname = "%s_ngseasy@%s.log" % (datetime.datetime.now().strftime("%d%m%y"), run_id)
    logfile = os.path.join(log_dir, logname)

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger("root")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(file_handler(logfile))
    logger.addHandler(console_handler())

    file_only = logging.getLogger("file-only")
    file_only.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_only.addHandler(file_handler(logfile))

    return logger


def multi_proc_handler(logfile):
    handler = MultiProcLogHandler(logfile)
    handler.setFormatter(FORMATTER)
    return handler


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

    logger = logging.getLogger("root")
    if len(logger.handlers) == 0:
        logger.addHandler(console_handler())
        logger.setLevel(logging.DEBUG)
    return logger


class MultiProcLogHandler(logging.Handler):
    def __init__(self, logfile):
        logging.Handler.__init__(self)

        self._handler = logging.FileHandler(logfile)
        self.queue = multiprocessing.Queue(-1)

        thrd = threading.Thread(target=self.receive)
        thrd.daemon = True
        thrd.start()

    def setFormatter(self, fmt):
        logging.Handler.setFormatter(self, fmt)
        self._handler.setFormatter(fmt)

    def receive(self):
        while True:
            try:
                record = self.queue.get()
                self._handler.emit(record)
            except (KeyboardInterrupt, SystemExit):
                raise
            except EOFError:
                break
            except:
                traceback.print_exc(file=sys.stderr)

    def send(self, s):
        self.queue.put_nowait(s)

    def _format_record(self, record):
        if record.args:
            record.msg = record.msg % record.args
            record.args = None
        if record.exc_info:
            dummy = self.format(record)
            record.exc_info = None

        return record

    def emit(self, record):
        try:
            s = self._format_record(record)
            self.send(s)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

    def close(self):
        self._handler.close()
        logging.Handler.close(self)
