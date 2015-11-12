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

import argparse

from logger import logger

import os.path


def parse_job_args(argv, name=""):
    uid = os.getuid()
    gid = os.getgid()

    parser = argparse.ArgumentParser(description="NGSPyEasy %s Job" % name)
    parser.add_argument("-c", "--config", dest="config", required=True, type=path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=existed_directory,
                        help="ngs_projects directory path")
    parser.add_argument("-i", "--sample_id", dest="sample_id", help="sample_id to run %s on" % name)
    parser.add_argument("-t", "--task", dest="task", default="no-task", help="%s task")
    parser.add_argument("-u", "--uid", dest="uid", type=int, default=uid, help="files owner uid")
    parser.add_argument("-g", "--gid", dest="gid", type=int, default=gid, help="files owner gid")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")
    args = parser.parse_args(argv)
    logger().debug("Parsed command line arguments: %s " % args)
    return args


def existed_directory(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError('%s is not an existed directory path' % path)
    return path


def existed_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError('%s is not an existed file path' % path)
    return path


def path_basename(path):
    return os.path.basename(path)
