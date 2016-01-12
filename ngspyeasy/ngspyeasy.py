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

import executor
import os
import cmdargs
import tsv_config
from logger import logger, init_main_logger
import yaml


def main(argv):
    parser = argparse.ArgumentParser(description="NGSpyeasy pipeline runner")
    parser.add_argument("playbook_path", metavar='/path/to/your_pipeline.yml', type=cmdargs.existed_file)
    parser.add_argument("--version", action="version", version="%(prog)s 3.0", help="print software version")
    parser.add_argument("--samples", metavar="/path/to/config.tsv", dest="samples_tsv", required=False,
                        type=cmdargs.existed_file, help="List of samples in TSV format")
    parser.add_argument("--vars", dest="var_files", metavar="/path/to/your/vars.yml", help="additional variables",
                        type=cmdargs.existed_file, action="append")
    parser.add_argument("--provider", dest="provider", choices=["local", "lsf"], default="local",
                        help="job scheduler")
    parser.add_argument("--log_dir", dest="log_dir", type=cmdargs.existed_directory)

    args = parser.parse_args(argv)

    if args.log_dir is not None:
        init_main_logger(args.log_dir)

    logger().debug("Command line arguments: %s" % args)

    samples_tsv = os.path.abspath(args.samples_tsv) if args.samples_tsv else None
    playbook_path = os.path.abspath(args.playbook_path)
    var_files = [os.path.abspath(f) for f in args.var_files]

    all_samples = read_samples(samples_tsv)
    logger().info("Number of samples: %s" % len(all_samples))

    all_tasks = read_tasks(playbook_path)
    logger().info("Number of tasks: %s" % len(all_tasks))

    vars = read_variables(var_files)
    vars["all_samples"] = all_samples

    options = []
    if args.log_dir is not None:
        options.append("--log_dir %s" % args.log_dir)
    if args.samples_tsv is not None:
        options.append("--samples %s" % args.samples_tsv)

    # TODO create structure to create sets of parallel tasks

    logger().info("Starting job scheduler: provider=%s" % args.provider)
    executor.start(provider=args.provider, log_dir=args.log_dir)

    try:
        for task_index, task in enumerate(all_tasks, start=0):
            jobs = []
            for name, cmd in task.as_commands(task_index, task, vars):
                executor.submit(name, cmd)
                jobs.append(name)
            while len(jobs) > 0:
                name = executor.results_queue.get()
                jobs.remove(name)

    except Exception as e:
        logger().exception(e)
    finally:
        executor.stop()


def read_tasks(playbook_path):
    logger().info("Reading playbook yaml...")
    with open(playbook_path, 'r') as stream:
        return yaml.load(stream)


def read_samples(samples_tsv):
    if samples_tsv is None:
        return []

    logger().debug("TSV config path: %s" % samples_tsv)
    tsv_conf = tsv_config.parse(samples_tsv)

    logger().info("TSV config first line: %s" % str(tsv_conf.row_at(0)))
    return list(tsv_conf.all_rows())


def read_variables(var_files):
    d = dict()
    for var_file in var_files:
        with open(var_file, 'r') as stream:
            vars = yaml.load(stream)
            d.update(vars)
    return d


def signal_handler(signum, frame):
    logger().info("Got SIGINT(%s) signal" % str(signum))
    executor.stop()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
