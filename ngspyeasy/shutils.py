#!/usr/bin/env python
import tempfile
from string import Template
import subprocess
import sys

import os
import stat


def run_command(cmd, logger):
    logger.debug("cmd to run: %s" % " ".join(cmd))
    proc = subprocess.Popen(
        ["/bin/bash", "-c",
         # ~/.bashrc can contain environment variables valuable for running the command; unfortunately it doesn't
         # run automatically in the docker container if bash runs in non-interactive mode (without '-i' flag). Read
         # the .bashrc source for more details.
         "python /ngspyeasy/bin/fix_bashrc.py; source ~/.bashrc_fixed; echo $CLASSPATH; " + " ".join(cmd)],
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
        logger.info("KeyboardInterrupt received")

    logger.debug("cmd:\n %s" % ''.join(lines))

    if proc.returncode:
        logger.error("Command [[\n%s\n]] failed. See logs for details", " ".join(cmd))


def script_from_template(template_path):
    if not os.path.isfile(template_path):
        raise IOError("Template file not found: %s" % template_path)

    with open(template_path) as f:
        lines = f.readlines()

    return ShellScript(lines)


class ShellScript(object):
    def __init__(self, lines):
        self.lines = [line.rstrip('\n') for line in lines]
        self.vars = dict()

    def add_variables(self, **kwargs):
        for key, value in kwargs.iteritems():
            self.vars[key] = value

    def variable_assignments(self):
        return ["%s=\"%s\"" % (key, value) for (key, value) in self.vars]

    def source(self):
        return self.lines[0:]

    def to_temporary_file(self, validate=True):
        if validate:
            t = Template("".join(self.lines))
            t.substitute(self.vars)

        source = ["#!/usr/bin/env bash"]
        source += self.variable_assignments()
        source += self.lines[1:] if len(self.lines) > 0 and self.lines[0].startswith("#!") else self.lines[0:]

        file = tempfile.NamedTemporaryFile()
        file.write("\n".join(source))
        file.close()

        st = os.stat(file.name)
        os.chmod(file.name, st.st_mode | stat.S_IEXEC)
        return file.name