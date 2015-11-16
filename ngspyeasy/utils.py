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

from functools import partial

def uniq_set(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


class Enum(object):
    def __init__(self, *keys):
        self.__dict__.update(zip(keys, range(len(keys))))


class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class LazyDict(dict):
    """
    A lazy dictionary implementation which will try
    to evaluate all values on access and cache the
    result for later access.
    """

    def set_lazy(self, key, item, *args, **kwargs):
        """
        Allow the setting of a callable and arguments
        as value of dictionary.
        """
        if callable(item):
            item = partial(item, *args, **kwargs)
        super(LazyDict, self).__setitem__(key, item)

    def __getitem__(self, key):
        item = super(LazyDict, self).__getitem__(key)
        try:
            self[key] = item = item()
        except TypeError:
            pass
        return item
