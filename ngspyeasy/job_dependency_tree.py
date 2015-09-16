#!/usr/bin/env python

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
        next_job = None
        queue = [self.root]
        visited = set()
        while queue:
            curr = queue.pop(0)
            visited.add(curr)

            if curr.is_new():
                next_job = curr
                break

            queue.extend(set(filter(lambda x: x.is_done() or x.is_new(), curr.get_children())) - visited)

        if next_job is not None:
            next_job.start()
            return next_job.get_id(), next_job.get_details()
        return None, None

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
