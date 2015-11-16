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
    def __init__(self, tmpl, tool_ref):
        self.name = tmpl.name
        self.template = tmpl.template
        self.vars = tmpl.vars
        self.inputs = tmpl.input
        self.outputs = tmpl.output
        self.tool_ref = tool_ref

    def files_must_exist(self, names, env):
        return self.files_exist(names, env, True)

    def files_exist(self, names, env, must=False):
        files = [env[x] for x in names]
        for f in files:
            if os.path.isfile(f):
                continue
            if must:
                raise ValueError("Can't find file: %s" % f)
            return False
        return True

    def key_values(self, env):
        d = dict()
        for k in self.vars + self.inputs + self.outputs:
            d[k] = env[k]
        return d

    def run(self, env):
        if self.files_exist(self.outputs, env):
            logger().info("Skipping this bit... results are already exist")
            return

        self.files_must_exist(self.inputs, env)

        tmpl = sh_template.load(self.tool_ref, self.template)
        shcmd.run_command(tmpl.create_sh_file(self.key_values(env)))


def find(p):
    tool_path = p.rsplit("#", 1)
    tool_path = [None] + tool_path if len(tool_path) == 1 else tool_path
    tool_name = tool_path[0]
    tool_ref = tool_path[1]

    root = os.path.dirname(__file__)
    tool_dir = os.path.join(root, "resources", tool_ref)
    logger().debug("Pipeline tool directory: %s" % tool_dir)

    if not os.path.isdir(tool_dir):
        return None

    main_json = os.path.join(tool_dir, "main.json")
    logger().debug("Path to main.json: %s" % main_json)

    with open(main_json, 'r') as stream:
        templates = json.load(stream)

    tmpl = templates[0] if tool_name is None else next(t for t in templates if t.name == tool_name)
    return PipelineTool(tmpl, tool_ref)
