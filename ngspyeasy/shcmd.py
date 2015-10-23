#!/usr/bin/env python
import functools
import subprocess
import sys
import itertools
from logger import logger

import os

def has_chmod_permissions(curr_uid, path):
    return curr_uid == 0 or curr_uid == os.stat(path).st_uid


def chmod(dir, dmode, fmode):
    has_permissions = functools.partial(has_chmod_permissions, os.getuid())
    for root, dirs, files in os.walk(dir):
        z = zip(dirs, itertools.repeat(dmode)) + zip(files, itertools.repeat(fmode))
        for p, mode in z:
            path = os.path.join(root, p)
            if has_permissions(path):
                os.chmod(path, mode)


def run_command(cmd):
    cmd2run = " ".join(cmd) if isinstance(cmd, list) else cmd
    logger().debug("cmd to run: %s" % cmd2run)
    proc = subprocess.Popen(
        ["/bin/bash", "-c",
         # ~/.bashrc can contain environment variables valuable for running the command; unfortunately it doesn't
         # run automatically in the docker container if bash runs in non-interactive mode (without '-i' flag). Read
         # the .bashrc source for more details.
         "python /ngspyeasy/bin/fix_bashrc.py; source ~/.bashrc_fixed; echo $CLASSPATH; " + cmd2run],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)

    lines = []
    try:
        for line in iter(proc.stdout.readline, b''):
            sys.stdout.write(line)
            sys.stdout.flush()
            lines.append(line)
        proc.stdout.close()
    except KeyboardInterrupt:
        logger().info("KeyboardInterrupt received")

    logger().debug("cmd:\n %s" % ''.join(lines))

    if proc.returncode > 0:
        raise IOError("Command [[\n%s\n]] failed. See logs for details" % cmd2run)
