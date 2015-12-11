#!/usr/bin/env python

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

def main():
    module = AnsibleModule(
        argument_spec=dict(
            command=dict(required=True, default=None, type='str'),
            image=dict(required=True, default=None, type='str'),
            volumes=dict(required=False, default=[], type='list'),
            environment=dict(required=False, default=[], type='list'),
            secure=dict(required=False, default=True, type='bool'),
            sudo=dict(required=False, default=True, type='bool'),
            rm=dict(required=False, default=True, type='bool')
        )
    )

    command = module.params['command']
    image = module.params['image']
    volumes = module.params['volumes']
    environment = module.params['environment']
    secure = module.params['secure']
    sudo = module.params['sudo']
    rm = module.params['rm']

    cmd = []
    if secure:
        cmd.append("sudo dockercmd run")
    else:
        cmd.append(("sudo " if sudo else "") + "docker run")

    if rm:
        cmd.append("--rm")

    cmd += ["-e " + x for x in environment]
    cmd += ["-v " + x for x in volumes]
    cmd += [image, command]

    proc = subprocess.Popen(
        ["/bin/bash", "-c", " ".join(cmd)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)

    lines = []

    for line in iter(proc.stdout.readline, b''):
        lines.append(line)
    proc.stdout.close()

    result = dict(stdout="".join(lines),
                  stdout_lines=lines)
    module.exit_json(**result)


# ===========================================

# import module snippets
from ansible.module_utils.basic import *

main()
