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
import pipeline_tools
import tsv_config
import projects_dir


def main(argv):
    parser = argparse.ArgumentParser(description="NGSPyEasy Tool")
    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory,
                        help="ngs_projects directory path")
    parser.add_argument("-r", "--resources-dir", dest="resources_dir", type=cmdargs.existed_directory,
                        help="ngs_resources directory path")
    parser.add_argument("-t", "--tool", dest="tool", required=True, help="pipeline tool name")
    parser.add_argument("-i", "--sample_id", dest="sample_id", help="sample_id to run with")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--version", action="version", version="%(prog)s 2.0", help="print software version")
    args = parser.parse_args(argv)
    logger().debug("Parsed command line arguments: %s " % args)

    projects_home = projects_dir.ProjectsDir(args.projects_dir, args.resources_dir)
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

    for row in rows2run(tsv_conf, args.sample_id):
        run_2l(args.tool, row, projects_home)


def rows2run(tsv_conf, sample_id):
    if sample_id is None:
        return tsv_conf.all_rows()
    return [x for x in tsv_conf.all_rows() if x.sample_id() == sample_id]


def run_2l(tool_path, row, projects_home):
    tool = pipeline_tools.find_tool(tool_path)
    if tool is None:
        raise ValueError("Unknown pipeline tool: %s" % tool)

    try:
        tool.run(row, projects_home)
    except Exception as ex:
        logger().exception(ex)


if __name__ == '__main__':
    main(sys.argv[1:])
