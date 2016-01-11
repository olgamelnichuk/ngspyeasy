#!/usr/bin/env python
import subprocess

import os


class Provider(object):
    def submit(self, name, cmd):
        raise NotImplementedError()

    def list(self):
        raise NotImplementedError()

    def terminate(self):
        raise NotImplementedError()


class LocalProvider(Provider):
    def __init__(self):
        self._processes = []

    def submit(self, name, cmd):
        proc = subprocess.Popen(["/bin/bash", "-c", cmd], env=os.environ.copy())
        self._processes.append((proc, cmd, name))
        return name

    def list(self):
        running = []
        for (proc, cmd, name) in self._processes:
            ret = proc.poll()
            if ret is None:
                running.append((proc, cmd, name))

        self._processes = running
        return [x[2] for x in running]

    def terminate(self, sigkill=False):
        unfinished = []
        for (proc, cmd, name) in self._processes:
            ret = proc.poll()
            if ret is None:
                unfinished.append((proc, cmd, name))
                if sigkill:
                    proc.kill()
                else:
                    proc.terminate()
        self._processes = unfinished


class LSFProvider(Provider):
    def __init__(self):
        pass
        # TODO
