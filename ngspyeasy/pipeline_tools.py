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
import os
import sh_template
import docker_env
import pipeline_env


class ToolSpec(object):
    def __init__(self, spec, tool_dir):
        self.spec = spec
        self.tool_dir = tool_dir

    def vars(self):
        return self.spec["vars"]

    def input(self):
        return self.spec["input"]

    def output(self):
        return self.spec["output"]

    def template(self):
        return self.spec["template"]

    def image(self):
        return self.spec["image"]

    def name(self):
        return self.spec["name"]

    def resource_path(self):
        return os.path.join(self.tool_dir, self.template())

    def id(self):
        return self.tool_dir + "#" + self.name()

    def all_vars(self):
        return self.vars() + self.input() + self.output()


class PipelineTool(object):
    def __init__(self, spec):
        self.spec = spec

    def output_files(self, env):
        return [env[x] for x in self.spec.output()]

    def input_files(self, env):
        return [env[x] for x in self.spec.input()]

    def files_must_exist(self, files):
        return self.files_exist(files, must=True)

    def files_exist(self, files, must=False):
        for f in files:
            if os.path.isfile(f):
                continue
            if must:
                raise ValueError("Can't find file: %s" % f)
            return False
        return True

    def all_vars(self, env):
        d = dict()
        for k in self.spec.all_vars():
            d[k] = env[k]
        return d

    def run(self, row, projects_home):
        host_env = pipeline_env.as_dict(row, projects_home)

        if self.files_exist(self.output_files(host_env)):
            logger().info("Skipping this bit... results are already exist")
            return

        self.files_must_exist(self.input_files(host_env))

        container_env = pipeline_env.as_dict(row, docker_env.projects_home())
        docker_env.run_command(self.cmd(container_env), self.spec.image(), projects_home)

    def cmd(self, env):
        tmpl = sh_template.load(self.spec.resource_path())
        return tmpl.as_executable(self.all_vars(env))


def normalize(p):
    arr = p.rsplit("#", 1)
    arr = arr + [None] if len(arr) == 1 else arr
    return arr[1], arr[0]


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
        specs = json.load(stream)

    return [ToolSpec(s, tool_dir) for s in specs]


def find_tool(p):
    tool_name, tool_dir = normalize(p)
    specs = find(p)

    if len(specs) == 0:
        return None

    spec = specs[0] if tool_name is None else next(s for s in specs if s.name == tool_name)
    return PipelineTool(ToolSpec(spec, tool_dir))
