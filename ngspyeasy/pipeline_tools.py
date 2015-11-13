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

import json

from logger import logger
import shcmd
import os
import sh_template


class PipelineTool(object):
    def __init__(self):
        # TODO
        pass

    def results_exist(self):
        # TODO
        return False

    def verify_inputs(self):
        # TODO
        pass

    def run(self, vars):
        if self.results_exist(vars):
            logger().info("Skipping this bit... results already exist")
            return

        self.verify_inputs(vars)

        tmpl = sh_template.load(path)
        shcmd.run_command(tmpl.create_sh_file(vars))


def find(p):
    root = os.path.dirname(__file__)
    tool_dir = os.path.join(root, "resources", p)
    logger().debug("Tool dir: %s" % tool_dir)

    if not os.path.isdir(tool_dir):
        return None

    main_json = os.path.join(tool_dir, "main.json")
    logger().debug("Path to main.json: %s" % main_json)

    with open(main_json, 'r') as stream:
        descr = json.load(stream)
