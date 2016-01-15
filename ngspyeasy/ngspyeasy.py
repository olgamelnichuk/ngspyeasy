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
import playbook_yaml
import os
import cmdargs
from logger import logger, init_main_logger


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

    playbook_path = os.path.abspath(args.playbook_path)
    samples_tsv = os.path.abspath(args.samples_tsv) if args.samples_tsv else None
    var_files = [os.path.abspath(f) for f in args.var_files]

    pb = playbook_yaml.parse(playbook_path, samples_tsv, var_files, args.log_dir)

    logger().info("Starting job executor: provider=%s" % args.provider)
    executor.start(provider=args.provider, log_dir=args.log_dir)

    try:
        for play in pb.plays():
            jobs = []
            for name, cmd in play.commands():
                executor.submit(name, cmd)
                jobs.append(name)
            while len(jobs) > 0:
                name = executor.results_queue.get()
                jobs.remove(name)

    except Exception as e:
        logger().exception(e)
    finally:
        executor.stop()


def signal_handler(signum, frame):
    logger().info("Got SIGINT(%s) signal" % str(signum))
    executor.stop()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)
    main(sys.argv[1:])
