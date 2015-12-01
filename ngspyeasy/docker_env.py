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
import sys

from logger import logger
import os
from docker import Client
import projects_dir

HOME = "/home/pipeman"

NGS_PROJECTS = HOME + "/ngs_projects"

NGS_RESOURCES = NGS_PROJECTS + "/ngseasy_resources"

DOCKER_OPTS = ""

DOCKER_BASEURL = 'unix://var/run/docker.sock'


def projects_home():
    return projects_dir.ProjectsDir(NGS_PROJECTS)


def volumes(projects_home):
    ngs_projects = projects_home.root()
    ngs_resources = projects_home.resources_dir()
    return [ngs_projects + ":" + NGS_PROJECTS,
            ngs_resources + ":" + NGS_RESOURCES]


def environment():
    return dict(HOME=HOME)


def working_dir():
    return HOME


def user():
    return str(os.getuid()) + ":" + str(os.getgid())


def run_command(cmd, image, name, projects_home):
    c = Client(base_url=DOCKER_BASEURL)
    container = c.create_container(image, command=cmd, name=name, volumes=volumes(projects_home),
                                   environment=environment(), working_dir=working_dir(), user=user())
    c.start(container, publish_all_ports=True)
    lines = []
    for line in c.logs(container, stdout=True, stderr=True, stream=True, timestamps=True, tail='all'):
        sys.stdout.write(line)
        sys.stdout.flush()
        lines.append(line)

    logger().debug("\n" + "".join(lines))

    status = c.wait(container)
    logger().debug("exit status code: %s" % status)
    c.remove_container(id)
    return status
