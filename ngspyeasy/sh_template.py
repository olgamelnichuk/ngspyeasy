#!/usr/bin/env python
import stat
from string import Template
import tempfile
from logger import log_debug
import os


def load(*p):
    root = os.path.dirname(__file__)
    path = os.path.join(root, *(["resources"] + list(p)))
    log_debug("Loading sh template from: %s" % path)

    if not os.path.isfile(path):
        raise IOError("Sh template file not found: %s" % path)

    with open(path) as f:
        lines = f.readlines()

    return ShTemplate(lines)


class ShTemplate(object):
    def __init__(self, lines):
        self.lines = [x for x in lines if not x.startswith("#")]

    def validate(self, **kwargs):
        t = Template("".join(self.lines))
        t.substitute(kwargs)

    def source(self, **kwargs):
        source = ["#!/usr/bin/env bash"]
        source += ["%s=\"%s\"" % (key, value) for (key, value) in kwargs.iteritems()]
        source += [Template("".join(self.lines)).safe_substitute()]
        return source

    def create_sh_file(self, validate=True, **kwargs):
        if validate:
            self.validate(**kwargs)

        source = self.source(**kwargs)
        log_debug("\n".join(source))

        file = tempfile.NamedTemporaryFile(delete=False)
        file.write("\n".join(source))
        file.close()

        st = os.stat(file.name)
        os.chmod(file.name, st.st_mode | stat.S_IEXEC)
        return file.name
