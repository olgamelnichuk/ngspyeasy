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

import os.path
import docker_env


def path_of_fixed_bashrc():
    bashrc = os.path.join(docker_env.HOME, ".bashrc")
    bashrc_fixed = os.path.join(docker_env.HOME, ".bashrc_fixed")
    if not os.path.isfile(bashrc):
        return bashrc_fixed

    with open(bashrc, 'r') as input:
        with open(bashrc_fixed, 'a') as output:
            lines = input.readlines()
            ignore = False
            for line in lines:
                output.write('#' + line if ignore else line)
                if line.startswith('# If not running interactively, don\'t do anything'):
                    ignore = True

                if line.startswith('esac') and ignore:
                    ignore = False
    return bashrc_fixed


if __name__ == '__main__':
    print path_of_fixed_bashrc()
