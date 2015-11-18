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


class ToolTemplate(object):
    def __init__(self, tmpl, tmpl_dir):
        self.tmpl = tmpl
        self.tmpl_dir = tmpl_dir

    def vars(self):
        return self.tmpl["vars"]

    def input(self):
        return self.tmpl["input"]

    def output(self):
        return self.tmpl["output"]

    def template(self):
        return self.tmpl["template"]

    def image(self):
        return self.tmpl["image"]

    def name(self):
        return self.tmpl["name"]

    def dir(self):
        return self.tmpl_dir

    def id(self):
        return self.dir() + "#" + self.name()

    def allvars(self):
        return self.vars() + self.input() + self.output()


class PipelineTool(object):
    def __init__(self, tmpl):
        self.tmpl = tmpl

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
        for k in self.tmpl.allvars():
            d[k] = env[k]
        return d

    def run(self, env):
        if self.files_exist(self.tmpl.output(), env):
            logger().info("Skipping this bit... results are already exist")
            return

        self.files_must_exist(self.tmpl.input(), env)

        tmpl = sh_template.load(self.tmpl.dir(), self.tmpl.template())
        shcmd.run_command(tmpl.create_sh_file(self.key_values(env)))


def normalize(p):
    tool_path = p.rsplit("#", 1)
    tool_path = [None] + tool_path if len(tool_path) == 1 else tool_path
    return tool_path[0], tool_path[1]


def tool_dirs():
    root = os.path.dirname(__file__)
    tools_root = os.path.join(root, "resources", "tools")

    for root, dirs, files in os.walk(tools_root):
        for f in files:
            if f == "main.json":
                yield os.path.relpath(root, tools_root)


def find(p):
    tool_name, tool_dir = normalize(p)

    root = os.path.dirname(__file__)
    absolute_dir = os.path.join(root, "resources", "tools", tool_dir)
    logger().debug("Pipeline tool directory: %s" % absolute_dir)

    if not os.path.isdir(absolute_dir):
        return []

    main_json = os.path.join(absolute_dir, "main.json")
    logger().debug("Path to main.json: %s" % main_json)

    with open(main_json, 'r') as stream:
        templates = json.load(stream)

    return [ToolTemplate(t, tool_dir) for t in templates]


def find_template(p):
    tool_name, tool_dir = normalize(p)
    templates = find(p)

    if len(templates) == 0:
        return None

    tmpl = templates[0] if tool_name is None else next(t for t in templates if t.name == tool_name)
    return PipelineTool(ToolTemplate(tmpl, tool_dir))
