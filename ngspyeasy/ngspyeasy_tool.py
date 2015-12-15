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
from ansible.module_utils import basic

import os
import cmdargs
from logger import logger, init_sample_logger
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
    parser.add_argument("--log_dir", dest="log_dir", type=cmdargs.existed_directory)

    args = parser.parse_args(argv)

    run_index = args.run_index
    run_id = args.run_id
    tsv_path = os.path.abspath(args.tsv_path)
    pipeline_path = os.path.abspath(args.pipeline_path)
    var_files = [os.path.abspath(f) for f in args.var_files]

    if args.log_dir is not None:
        init_sample_logger(args.log_dir, run_id)

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

        library_dir = os.path.join(os.path.dirname(pipeline_path), "library")
        if os.path.exists(library_dir):
            shutil.copytree(library_dir, os.path.join(temp_dir, "library"))

        playbook = os.path.join(temp_dir, "playbook.yml")
        with open(playbook, 'w') as outfile:
            outfile.write(yaml.dump([task], default_flow_style=False))

        run_playbook(temp_dir, extra_vars)
    finally:
        shutil.rmtree(temp_dir)


def run_playbook(dir, extra_vars):
    utils.VERBOSITY = 0
    playbook_cb = MyPlaybookCallbacks(verbose=utils.VERBOSITY)
    stats = callbacks.AggregateStats()
    runner_cb = MyPlaybookRunnerCallbacks(stats, verbose=utils.VERBOSITY)

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
    logger().info(results)


def file_logger():
    return logger(file_only=True)


class MyPlaybookCallbacks(callbacks.PlaybookCallbacks):
    def __init__(self, verbose=False):
        super(MyPlaybookCallbacks, self).__init__(verbose)

    def on_start(self):
        super(MyPlaybookCallbacks, self).on_start()

    def on_notify(self, host, handler):
        super(MyPlaybookCallbacks, self).on_notify(host, handler)

    def on_no_hosts_matched(self):
        super(MyPlaybookCallbacks, self).on_no_hosts_matched()

    def on_no_hosts_remaining(self):
        super(MyPlaybookCallbacks, self).on_no_hosts_remaining()

    def on_task_start(self, name, is_conditional):
        file_logger().info("task_start: %s" % name)
        super(MyPlaybookCallbacks, self).on_task_start(name, is_conditional)

    def on_vars_prompt(self, varname, private=True, prompt=None, encrypt=None, confirm=False, salt_size=None, salt=None,
                       default=None):
        return super(MyPlaybookCallbacks, self).on_vars_prompt(varname, private, prompt, encrypt, confirm, salt_size,
                                                               salt, default)

    def on_setup(self):
        file_logger().info("GATHERING FACTS...")
        super(MyPlaybookCallbacks, self).on_setup()

    def on_import_for_host(self, host, imported_file):
        file_logger().info("%s: importing file: %s" % (host, imported_file))
        super(MyPlaybookCallbacks, self).on_import_for_host(host, imported_file)

    def on_not_import_for_host(self, host, missing_file):
        file_logger().info("%s: not importing file: %s" % (host, missing_file))
        super(MyPlaybookCallbacks, self).on_not_import_for_host(host, missing_file)

    def on_play_start(self, name):
        file_logger().info("PLAY [%s]" % name)
        super(MyPlaybookCallbacks, self).on_play_start(name)

    def on_stats(self, stats):
        super(MyPlaybookCallbacks, self).on_stats(stats)


class MyPlaybookRunnerCallbacks(callbacks.PlaybookRunnerCallbacks):
    def __init__(self, stats, verbose=None):
        super(MyPlaybookRunnerCallbacks, self).__init__(stats, verbose)

    def on_unreachable(self, host, results):
        if self.runner.delegate_to:
            host = '%s -> %s' % (host, self.runner.delegate_to)

        item = None
        if type(results) == dict:
            item = results.get('item', None)
            if isinstance(item, unicode):
                item = utils.unicode.to_bytes(item)
            results = basic.json_dict_unicode_to_bytes(results)
        else:
            results = utils.unicode.to_bytes(results)
        host = utils.unicode.to_bytes(host)
        if item:
            msg = "fatal: [%s] => (item=%s) => %s" % (host, item, results)
        else:
            msg = "fatal: [%s] => %s" % (host, results)
        file_logger().error(msg)
        super(MyPlaybookRunnerCallbacks, self).on_unreachable(host, results)

    def on_failed(self, host, results, ignore_errors=False):
        if self.runner.delegate_to:
            host = '%s -> %s' % (host, self.runner.delegate_to)

        results2 = results.copy()
        results2.pop('invocation', None)

        item = results2.get('item', None)
        parsed = results2.get('parsed', True)
        returned_msg = results2.pop('msg', None)
        module_msg = ''
        if not parsed:
            module_msg = results2.pop('msg', None)

        if item:
            msg = "failed: [%s] => (item=%s) => %s" % (host, item, utils.jsonify(results2))
        else:
            msg = "failed: [%s] => %s" % (host, utils.jsonify(results2))
        file_logger().error(msg)

        if returned_msg:
            file_logger().error(returned_msg)
        if not parsed and module_msg:
            file_logger().error(module_msg)
        if ignore_errors:
            file_logger().info("...ignoring")

        super(MyPlaybookRunnerCallbacks, self).on_failed(host, results, ignore_errors=ignore_errors)

    def on_ok(self, host, host_result):
        item = host_result.get('item', None)

        host_result2 = host_result.copy()
        host_result2.pop('invocation', None)
        verbose_always = host_result2.pop('verbose_always', False)
        changed = host_result.get('changed', False)
        ok_or_changed = 'ok'
        if changed:
            ok_or_changed = 'changed'

        # show verbose output for non-setup module results if --verbose is used
        msg = ''
        if (not self.verbose or host_result2.get("verbose_override", None) is not
            None) and not verbose_always:
            if item:
                msg = "%s: [%s] => (item=%s)" % (ok_or_changed, host, item)
            else:
                if 'ansible_job_id' not in host_result or 'finished' in host_result:
                    msg = "%s: [%s]" % (ok_or_changed, host)
        else:
            # verbose ...
            if item:
                msg = "%s: [%s] => (item=%s) => %s" % (
                ok_or_changed, host, item, utils.jsonify(host_result2, format=verbose_always))
            else:
                if 'ansible_job_id' not in host_result or 'finished' in host_result2:
                    msg = "%s: [%s] => %s" % (ok_or_changed, host, utils.jsonify(host_result2, format=verbose_always))

        if msg != '':
            file_logger().info(msg)
        if 'warnings' in host_result2 and host_result2['warnings']:
            for warning in host_result2['warnings']:
                file_logger().warn("warning: %s" % warning)
        super(MyPlaybookRunnerCallbacks, self).on_ok(host, host_result)

    def on_skipped(self, host, item=None):
        if self.runner.delegate_to:
            host = '%s -> %s' % (host, self.runner.delegate_to)
        if item:
            msg = "skipping: [%s] => (item=%s)" % (host, item)
        else:
            msg = "skipping: [%s]" % host
        file_logger().info(msg)
        super(MyPlaybookRunnerCallbacks, self).on_skipped(host, item)

    def on_no_hosts(self):
        file_logger().error("FATAL: no hosts matched or all hosts have already failed -- aborting\n")
        super(MyPlaybookRunnerCallbacks, self).on_no_hosts()

    def on_async_poll(self, host, res, jid, clock):
        super(MyPlaybookRunnerCallbacks, self).on_async_poll(host, res, jid, clock)

    def on_async_ok(self, host, res, jid):
        if jid:
            msg = "<job %s> finished on %s" % (jid, host)
            file_logger().info(msg)
        super(MyPlaybookRunnerCallbacks, self).on_async_ok(host, res, jid)

    def on_async_failed(self, host, res, jid):
        msg = "<job %s> FAILED on %s" % (jid, host)
        file_logger().error(msg)
        super(MyPlaybookRunnerCallbacks, self).on_async_failed(host, res, jid)

    def on_file_diff(self, host, diff):
        file_logger().info(utils.get_diff(diff))
        super(MyPlaybookRunnerCallbacks, self).on_file_diff(host, diff)


if __name__ == '__main__':
    main(sys.argv[1:])
