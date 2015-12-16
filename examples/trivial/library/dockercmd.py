#!/usr/bin/env python

def main():
    module = AnsibleModule(
        argument_spec=dict(
            command=dict(required=True, default=None, type='str'),
            image=dict(required=True, default=None, type='str'),
            volumes=dict(required=False, default=[], type='list'),
            environment=dict(required=False, default=[], type='list'),
            working_dir=dict(required=False, default=None, type='str'),
            secure=dict(required=False, default=True, type='bool'),
            sudo=dict(required=False, default=True, type='bool'),
            rm=dict(required=False, default=True, type='bool')
        )
    )

    command = module.params['command']
    image = module.params['image']
    volumes = module.params['volumes']
    environment = module.params['environment']
    working_dir = module.params['working_dir']
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

    if working_dir is not None:
        cmd.append("-w " + working_dir)

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
