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

import os
import cmdargs
import job_id_generator
import job_scheduler
import tsv_config
from logger import logger, init_main_logger
import yaml
import jinja2

TEST_MODE = False


class JobSubmitter(object):
    def __init__(self, mode="local", log_dir=None):
        self.dependencies = dict()
        self.mode = mode
        self.log_dir = log_dir

    def dependencies_for(self, sample_index):
        return [x for x in [self.dependencies.get(str(sample_index), None)] if x is not None]

    def update_dependencies(self, sample_index, job_id):
        self.dependencies[str(sample_index)] = job_id

    def submit(self, sample_index, task_index, config_path, pipeline_script, var_files):
        job_cmd = self.cmd(sample_index, task_index, config_path, pipeline_script, var_files)

        job_id = job_id_generator.get_next([str(sample_index), str(task_index)])
        job_dependencies = self.dependencies_for(sample_index)

        logger().debug(
            "Submit job(job_id='%s', dependencies='%s')" % (job_id, job_dependencies))

        job_scheduler.submit(
            job_id, self.create_cmd(job_id, job_dependencies, job_cmd), job_dependencies)
        self.update_dependencies(sample_index, job_id)

    def create_cmd(self, job_id, job_dependencies, job_cmd):
        return self.lsf_wrap(job_id, job_dependencies, job_cmd) if self.in_lsf_mode() else job_cmd

    def lsf_wrap(self, job_id, job_dependencies, cmd):
        lsf_dep_expression = ""
        if len(job_dependencies) > 0:
            lsf_dep_expression = "-w \"%s\"" % " && ".join(["ended(%s)" % x for x in job_dependencies])

        return "bsub -J %s %s %s" % (job_id, lsf_dep_expression, cmd)

    def in_lsf_mode(self):
        return self.mode == "lsf"

    def cmd(self, sample_index, task_index, config_path, pipeline_script, var_files):
        executable = "ngspyeasy_tool"
        if TEST_MODE:
            root = os.path.dirname(__file__)
            executable = "python " + os.path.abspath(os.path.join(root, "ngspyeasy_tool.py"))
        log_dir = []
        if self.log_dir:
            log_dir = ["--log_dir", self.log_dir]
        return " ".join([executable,
                         pipeline_script,
                         "--sample_index", str(sample_index),
                         "--task_index", str(task_index),
                         "--samples", config_path
                         ] + ["--vars " + x for x in var_files] + log_dir)


def main(argv):
    parser = argparse.ArgumentParser(description="NGSpeasy pipelines")
    parser.add_argument("pipeline_path", metavar='/path/to/your_pipeline.yml', type=cmdargs.existed_file)
    parser.add_argument("--version", action="version", version="%(prog)s 3.0", help="print software version")
    parser.add_argument("--samples", metavar="/path/to/config.tsv", dest="tsv_path", required=True,
                        type=cmdargs.existed_file, help="List of samples in TSV format")
    parser.add_argument("--vars", dest="var_files", metavar="/path/to/your/vars.yml", help="additional variables",
                        type=cmdargs.existed_file, action="append")
    parser.add_argument("--mode", dest="mode", choices=["local", "lsf"], default="local", help="job scheduler")
    parser.add_argument("--log_dir", dest="log_dir", type=cmdargs.existed_directory)

    args = parser.parse_args(argv)

    if args.log_dir is not None:
        init_main_logger(args.log_dir)

    tsv_path = os.path.abspath(args.tsv_path)
    pipeline_path = os.path.abspath(args.pipeline_path)
    var_files = [os.path.abspath(f) for f in args.var_files]

    logger().debug("Command line arguments: %s" % args)
    logger().debug("TSV config path: %s" % tsv_path)
    try:
        tsv_conf = tsv_config.parse(tsv_path)
    except (IOError, ValueError) as e:
        logger().exception(e)
        sys.exit(1)

    logger().info("TSV config first line: %s" % str(tsv_conf.row_at(0)))
    all_samples = tsv_conf.all_rows()
    logger().info("Number of samples: %s" % len(all_samples))

    try:
        with open(pipeline_path, 'r') as stream:
            tasks = yaml.load(stream)
    except yaml.scanner.ScannerError, e:
        logger().exception(e)
        sys.exit(1)

    try:
        logger().info("Starting job scheduler...")
        scheduler = job_scheduler.JobScheduler()
        scheduler.start()
    except Exception, e:
        logger().exception(e)
        sys.exit(1)

    try:
        submitter = JobSubmitter(mode=args.mode, log_dir=args.log_dir)
        for task_index, task in enumerate(tasks, start=0):
            tmpl = task["samples"]
            samples2run = all_samples
            if tmpl != "all":
                samples_str = jinja2.Template(tmpl).render({"all_samples": all_samples})
                samples2run = eval(samples_str)

            for sample_index, sample in enumerate(samples2run, start=0):
                submitter.submit(sample_index, task_index, tsv_path, pipeline_path, var_files)

        job_scheduler.all_done()
    except Exception as e:
        logger().exception(e)
        job_scheduler.stop()

    while True:
        threads = threading.enumerate()
        if len(threads) == 1:
            break
        for t in threads:
            if t != threading.currentThread():
                t.join(1)


def signal_handler(signum, frame):
    logger().info("Got SIGINT(%s) signal" % str(signum))
    job_scheduler.stop()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
