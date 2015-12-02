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

import stat
from string import Template
import tempfile
from logger import logger
import os


def load(*p):
    root = os.path.dirname(__file__)
    path = os.path.join(root, *(["resources"] + list(p)))
    return ShTemplate(path)


class ShTemplate(object):
    def __init__(self, path):
        logger().debug("Loading sh template from: %s" % path)

        if not os.path.isfile(path):
            raise IOError("sh template file not found: %s" % path)

        with open(path) as f:
            lines = f.readlines()

        self.lines = [x for x in lines if not x.startswith("#")]
        self.path = path

    def validate(self, **kwargs):
        try:
            t = Template("".join(self.lines))
            t.substitute(kwargs)
        except KeyError, e:
            logger().exception(e)
            raise ValueError("SH Template validation failure (var: %s): %s" % (e.message, self.path))

    def source(self, **kwargs):
        source = ["#!/usr/bin/env bash"]
        source += ["%s=\"%s\"" % (key, value) for (key, value) in kwargs.iteritems()]
        source += [Template("".join(self.lines)).safe_substitute()]
        return source

    def as_executable(self, tmpdir, validate=True, **kwargs):
        if validate:
            self.validate(**kwargs)

        source = self.source(**kwargs)
        logger().debug("\n".join(source))

        file = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        file.write("\n".join(source))
        file.close()

        st = os.stat(file.name)
        os.chmod(file.name, st.st_mode | stat.S_IEXEC)
        return file.name
