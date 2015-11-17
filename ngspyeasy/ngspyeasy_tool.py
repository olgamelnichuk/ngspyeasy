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
import sys

import cmdargs
from logger import logger, init_logger
from ngspyeasy import pipeline_environment
import pipeline_tools
import tsv_config
import projects_dir
import os


def fix_file_permissions(projects_home, row, uid, gid):
    projects_home.fix_file_permissions(row.project_id(), row.sample_id(), uid, gid)


def main(argv):
    uid = os.getuid()
    gid = os.getgid()

    parser = argparse.ArgumentParser(description="NGSPyEasy Tool")
    parser.add_argument("-t", "--tool", dest="tool", required=True, help="relative tool path")
    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory,
                        help="ngs_projects directory path")
    parser.add_argument("-i", "--sample_id", dest="sample_id", help="sample_id to run with")
    parser.add_argument("-u", "--uid", dest="uid", type=int, default=uid, help="files owner uid")
    parser.add_argument("-g", "--gid", dest="gid", type=int, default=gid, help="files owner gid")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 2.0", help="print software version")
    args = parser.parse_args(argv)
    logger().debug("Parsed command line arguments: %s " % args)

    projects_home = projects_dir.ProjectsDir(args.projects_dir)
    log_file = projects_home.sample_log_file(args.config, args.sample_id)
    print "Opening log file: %s" % log_file

    init_logger(log_file, args.verbose)
    logger().info("Starting...")
    logger().debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    logger().info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        logger().error(e)
        sys.exit(1)

    try:
        for row in rows2run(tsv_conf, args.sample_id):
            run_2l(args.tool, row, projects_home, args.uid, args.gid)
    except Exception as ex:
        logger().exception(ex)
        sys.exit(1)


def rows2run(tsv_conf, sample_id):
    if sample_id is None:
        return tsv_conf.all_rows()
    return [x for x in tsv_conf.all_rows() if x.sample_id() == sample_id]


def run_2l(tool_path, row, projects_home, uid, gid):
    tool = pipeline_tools.find_template(tool_path)
    if tool is None:
        raise ValueError("Unknown pipeline tool: %s" % tool)

    try:
        tool.run(pipeline_environment.as_dict(row, projects_home))
    except Exception as ex:
        logger().exception(ex)
    finally:
        fix_file_permissions(projects_home, row, uid, gid)


if __name__ == '__main__':
    main(sys.argv[1:])
