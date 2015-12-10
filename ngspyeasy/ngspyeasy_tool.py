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
import shutil
import sys
import tempfile

import os
import cmdargs
from logger import logger
import tsv_config
from ansible.playbook import PlayBook
from ansible import callbacks
from ansible import utils
import yaml


def main(argv):
    parser = argparse.ArgumentParser(description="NGSpeasy pipelines")
    parser.add_argument("pipeline_path", metavar='/path/to/your_pipeline.yml', type=cmdargs.existed_file)
    parser.add_argument("--run_id", dest="run_id", help="sample_id to run with")
    parser.add_argument("--run_index", dest="run_index", type=int, help="task index to run")
    parser.add_argument("--version", action="version", version="%(prog)s 3.0", help="print software version")
    parser.add_argument("--samples", metavar="/path/to/config.tsv", dest="tsv_path", required=True,
                        type=cmdargs.existed_file, help="List of samples in TSV format")
    parser.add_argument("--vars", dest="var_files", metavar="/path/to/your/vars.yml", help="additional variables",
                        type=cmdargs.existed_file, action="append")

    args = parser.parse_args(argv)

    run_index = args.run_index
    run_id = args.run_id
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

    extra_vars = dict()

    sample = next(x for x in tsv_conf.all_rows() if x.sample_id == run_id)
    if sample is None:
        raise ValueError("Unknown sample id: %s" % run_id)

    extra_vars["sample"] = sample

    for var_file in var_files:
        with open(var_file, 'r') as stream:
            vars = yaml.load(stream)
            print vars
            extra_vars.update(vars)

    try:
        with open(pipeline_path, 'r') as stream:
            tasks = yaml.load(stream)
        task = tasks[run_index]
    except yaml.scanner.ScannerError, e:
        logger().exception(e)
        sys.exit(1)

    task.pop("samples", None)
    task["hosts"] = "all"

    temp_dir = tempfile.mkdtemp()
    print temp_dir
    try:
        roles_dir = os.path.join(os.path.dirname(pipeline_path), "roles")
        if os.path.exists(roles_dir):
            shutil.copytree(roles_dir, os.path.join(temp_dir, "roles"))

        playbook = os.path.join(temp_dir, "playbook.yml")
        with open(playbook, 'w') as outfile:
            outfile.write(yaml.dump([task], default_flow_style=False))

        run_playbook(temp_dir, extra_vars)
    finally:
        shutil.rmtree(temp_dir)


def run_playbook(dir, extra_vars):
    utils.VERBOSITY = 0
    playbook_cb = callbacks.PlaybookCallbacks(verbose=utils.VERBOSITY)
    stats = callbacks.AggregateStats()
    runner_cb = callbacks.PlaybookRunnerCallbacks(stats, verbose=utils.VERBOSITY)

    inventory = """
[localhost]
localhost ansible_connection=local
"""

    # Create a temporary file and write the template string to it
    hosts = tempfile.NamedTemporaryFile(delete=False, dir=dir)
    hosts.write(inventory)
    hosts.close()

    pb = PlayBook(
        playbook=os.path.join(dir, "playbook.yml"),
        host_list=hosts.name,
        callbacks=playbook_cb,
        runner_callbacks=runner_cb,
        extra_vars=extra_vars,
        stats=stats
    )

    results = pb.run()

    # Ensure on_stats callback is called
    # for callback modules
    playbook_cb.on_stats(pb.stats)
    print results


if __name__ == '__main__':
    main(sys.argv[1:])
