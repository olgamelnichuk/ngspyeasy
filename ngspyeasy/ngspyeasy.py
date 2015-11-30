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
import signal
import threading
import itertools

import os
import cmdargs
import job_id_generator
import pipeline_tools
import projects_dir
import job_scheduler
import tsv_config
from logger import init_logger, logger


class JobSubmitter(object):
    def __init__(self, projects_home, mode="docker"):
        self.dependencies = dict()
        self.projects_home = projects_home
        self.mode = mode

    def dependencies_for(self, sample_id):
        return [x for x in [self.dependencies.get(sample_id, None)] if x is not None]

    def update_dependencies(self, sample_id, job_id):
        self.dependencies[sample_id] = job_id

    def submit(self, tool, sample_id, config_name, verbose):
        job_cmd = self.cmd(tool, sample_id, config_name, verbose)
        tag = tool.name()

        job_id = job_id_generator.get_next([tag, sample_id])
        job_dependencies = self.dependencies_for(sample_id)

        logger().debug(
            "Submit job(job_id='%s', dependencies='%s')" % (job_id, job_dependencies))

        job_scheduler.submit(
            job_id, self.create_cmd(job_id, job_dependencies, job_cmd), job_dependencies)
        self.update_dependencies(sample_id, job_id)

    def create_cmd(self, job_id, job_dependencies, job_cmd):
        return self.lsf_wrap(job_id, job_dependencies, job_cmd) if self.in_lsf_mode() else job_cmd

    def lsf_wrap(self, job_id, job_dependencies, cmd):
        lsf_dep_expression = ""
        if len(job_dependencies) > 0:
            lsf_dep_expression = "\"%s\"" % " && ".join(["ended(%s)" % x for x in job_dependencies])

        return "bsub -J %s -w %s %s" % (job_id, lsf_dep_expression, cmd)

    def in_lsf_mode(self):
        return self.mode == "docker-lsf"

    def cmd(self, tool, sample_id, config_name, verbose):
        # root = os.path.dirname(__file__)
        # executable = os.path.join(root, "ngspyeasy_tool.py")
        return " ".join(["ngspyeasy_tool",
                         "-c", config_name,
                         "-d", self.projects_home.root(),
                         "-r", self.projects_home.resources_dir(),
                         "--sample_id", sample_id,
                         "--tool", tool.id(),
                         "--verbose" if verbose else ""])


def main(argv):
    parser = argparse.ArgumentParser(description="Python version of NGSeasy pipeline.")
    parser.add_argument("init", nargs="?")
    parser.add_argument("fastqc", nargs="?")
    parser.add_argument("trimmomatic", nargs="?")
    parser.add_argument("alignment", nargs="?")
    parser.add_argument("realign", nargs="?")
    parser.add_argument("bsqr", nargs="?")
    parser.add_argument("variant_calling", nargs="?")

    parser.add_argument("-c", "--config", dest="config", required=True, type=cmdargs.path_basename,
                        help="TSV configuration file name")
    parser.add_argument("-d", "--projects-dir", dest="projects_dir", required=True, type=cmdargs.existed_directory,
                        help="ngs_projects directory path")
    parser.add_argument("-r", "--resources-dir", dest="resources_dir", type=cmdargs.existed_directory,
                        help="ngs_resources directory path")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="turn ON verbose mode")
    parser.add_argument("--test", dest="test", action="store_true", help="turn ON test mode")
    parser.add_argument("-m", "--mode", dest="mode", choices=["docker", "docker-lsf"], default="docker",
                        help="running mode")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1", help="print software version")

    args = parser.parse_args(argv)

    projects_home = projects_dir.ProjectsDir(args.projects_dir, args.resources_dir)
    log_file = projects_home.main_log_file(args.config)
    print "Using log file: %s" % log_file

    logger = init_logger(log_file, args.verbose)
    logger.info("Starting NGSPyEasy 2.0...")
    logger.debug("Command line arguments: %s" % args)

    tsv_config_path = projects_home.config_path(args.config)
    logger.info("Reading TSV config: %s" % tsv_config_path)
    try:
        tsv_conf = tsv_config.parse(tsv_config_path)
    except (IOError, ValueError) as e:
        logger.error(e)
        sys.exit(1)

    logger.info("TSV config stats: %s" % tsv_conf.stats())

    if not args.init:
        try:
            logger.info("Starting job scheduler...")
            scheduler = job_scheduler.JobScheduler(args.test)
            scheduler.start()
        except Exception, e:
            logger.exception(e)
            sys.exit(1)

    verbose = args.verbose
    try:
        if init_required(args):
            ngspyeasy_init(tsv_conf, projects_home)

        submitter = JobSubmitter(projects_home, args.mode)
        for tool, sample_id in command_list(tsv_conf, args):
            submitter.submit(tool, sample_id, tsv_conf.filename(), verbose)
        job_scheduler.all_done()
    except Exception as e:
        logger.exception(e)
        job_scheduler.stop()

    while True:
        threads = threading.enumerate()
        if len(threads) == 1:
            break
        for t in threads:
            if t != threading.currentThread():
                t.join(1)


def init_required(args):
    if args.init:
        return True
    return not (
        args.fastqc and
        args.trimmomatic and
        args.alignment and
        args.realign and
        args.bsqr and
        args.variant_calling)


def ngspyeasy_init(tsv_conf, projects_home):
    logger().info("Checking ngs_projects directory structure...")
    projects_home.init_structure(tsv_conf)

    logger().info("Checking raw FastQ files...")
    projects_home.check_fastq(tsv_conf)


def command_list(tsv_conf, args):
    if args.init:
        return list()

    if args.fastqc:
        return ngspyeasy_fastqc(tsv_conf)

    if args.trimmomatic:
        return ngspyeasy_trimmomatic(tsv_conf)

    if args.alignment:
        return ngspyeasy_alignment(tsv_conf)

    if args.realign:
        return ngspyeasy_realn(tsv_conf)

    if args.bsqr:
        return ngspyeasy_bsqr(tsv_conf)

    if args.variant_calling:
        return ngspyeasy_variant_calling(tsv_conf)

    return ngspyeasy(tsv_conf)


def ngspyeasy(tsv_conf):
    return itertools.chain(ngspyeasy_fastqc(tsv_conf),
                           ngspyeasy_trimmomatic(tsv_conf),
                           ngspyeasy_alignment(tsv_conf),
                           ngspyeasy_realn(tsv_conf),
                           ngspyeasy_bsqr(tsv_conf),
                           ngspyeasy_variant_calling(tsv_conf))


def ngspyeasy_fastqc(tsv_conf):
    logger().info("Generating FastQC jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        fastqc_type = row.fastqc()

        if fastqc_type == "no-fastqc":
            logger().info("[%s] No fastqc jobs to be run for sample: '%s'" % (fastqc_type, sample_id))
            continue

        tools = pipeline_tools.find(os.path.join("fastqc", fastqc_type))

        if len(tools) == 0:
            raise ValueError("Unknown fastqc type: %s" % fastqc_type)

        for tmpl in tools:
            yield tmpl, sample_id


def ngspyeasy_trimmomatic(tsv_conf):
    logger().info("Generating Trimmomatic jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        trim_type = row.trim()

        if trim_type == "no-trim":
            logger().info(
                "[%s] No trimmomatic jobs to be run for sample: '%s'. NOT RECOMMENDED" % (trim_type, sample_id))
            continue

        tools = pipeline_tools.find(os.path.join("trimmomatic", trim_type))
        common = pipeline_tools.find(os.path.join("trimmomatic", "__after__"))

        if len(tools) == 0:
            raise ValueError("Unknown trimmomatic type: %s" % trim_type)

        for tool in (tools + common):
            yield tool, sample_id


def ngspyeasy_alignment(tsv_conf):
    logger().info("Generating Alignment jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        aligner_type = row.aligner()

        if aligner_type == "no-align":
            logger().info("[%s] No alignment jobs to be run for sample: '%s'." % (aligner_type, sample_id))
            continue

        tools = pipeline_tools.find(os.path.join("alignment", aligner_type))

        if len(tools) == 0:
            raise ValueError("Unknown aligner type: %s" % aligner_type)

        for tool in tools:
            yield tool, sample_id


def ngspyeasy_realn(tsv_conf):
    logger().info("Generating Realignment jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        realn_type = row.realn()

        if realn_type == "no-realn":
            logger().info("[%s] Skipping Indel Realignment for sample: '%s'." % (realn_type, sample_id))
            continue

        tools = pipeline_tools.find(os.path.join("realn", realn_type))

        if len(tools) == 0:
            raise ValueError("Unknown realign type: %s" % realn_type)

        for tool in tools:
            yield tool, sample_id


def ngspyeasy_bsqr(tsv_conf):
    logger().info("Generating Base Score Quality Recalibration jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        bsqr_type = row.bsqr()

        if bsqr_type == "no-bsqr":
            logger().info("[%s] Skipping Base quality score recalibration for sample: '%s'" % (bsqr_type, sample_id))
            continue

        tools = pipeline_tools.find(os.path.join("bsqr", bsqr_type))

        if len(tools) == 0:
            raise ValueError("Unknown bsqr type: %s" % bsqr_type)

        for tool in tools:
            yield tool, sample_id


def ngspyeasy_variant_calling(tsv_conf):
    logger().info("Generating Variant Calling jobs...")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        vc_type = row.varcaller()

        if vc_type == "no-vc":
            logger().info("[%s] Skipping Variant Calling for sample: '%s'" % (vc_type, sample_id))
            continue

        common = pipeline_tools.find(os.path.join("vc", "__before__"))
        tools = pipeline_tools.find(os.path.join("vc", vc_type))

        if len(tools) == 0:
            raise ValueError("Unknown variant calling type: %s" % vc_type)

        for tool in (common + tools):
            yield tool, sample_id


def signal_handler(signum, frame):
    logger().info("Got SIGINT(%s) signal" % str(signum))
    job_scheduler.stop()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
