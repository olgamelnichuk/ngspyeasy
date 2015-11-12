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

from utils import Enum

states = Enum('NEW', 'DONE', 'ERROR', 'RUNNING')


class JobDependencyTree(object):
    def __init__(self):
        self.root = Job("root", None)
        self.root.finish(0)

        self.index = 0
        self.dict = {"root": self.root}

    def append(self, job_id, job_dependencies=None, job_details=None):
        if self.dict.has_key(job_id):
            raise ValueError("Duplicated job_id: [%s]", job_id)

        if job_dependencies is None or len(job_dependencies) == 0:
            job_dependencies = ["root"]

        job = Job(job_id, job_details)
        parents = []
        for parent_id in job_dependencies:
            parent = self.dict.get(parent_id)
            if parent is None:
                raise ValueError("Job with id=[%s] doesn't exist", parent_id)
            parents.append(parent)

        self.dict[job_id] = job
        for parent in parents:
            parent.append_child(job)

    def get(self):
        next_job = self.bfs(lambda x: x.is_new())
        if next_job is not None:
            next_job.start()
            return next_job.get_id(), next_job.get_details()
        return None, None

    def has_running_jobs(self):
        job = self.bfs(lambda x: x.is_running())
        return job is not None

    def bfs(self, predicate):
        queue = [self.root]
        visited = set()
        while queue:
            curr = queue.pop(0)
            visited.add(curr)

            if predicate(curr):
                return curr

            queue.extend(set([x for x in curr.get_children() if x.is_done() or predicate(x)]) - visited)

        return None

    def done(self, job_id, retcode):
        job = self.dict.get(job_id)
        if job is None:
            raise ValueError("Job with id=[%s] doesn't exist", job_id)
        job.finish(retcode)

    def has_key(self, req_id):
        return self.dict.has_key(req_id)


class Job():
    def __init__(self, id, details):
        self.state = states.NEW
        self.children = set()
        self.id = id
        self.details = details

    def finish(self, retcode):
        self.state = states.DONE if retcode == 0 else states.ERROR

    def start(self):
        self.state = states.RUNNING

    def append_child(self, child_job):
        self.children.add(child_job)

    def get_id(self):
        return self.id

    def get_children(self):
        return set(self.children)

    def get_details(self):
        return self.details

    def is_new(self):
        return self.state == states.NEW

    def is_done(self):
        return self.state == states.DONE

    def is_running(self):
        return self.state == states.RUNNING

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.id == other.id
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.id)
