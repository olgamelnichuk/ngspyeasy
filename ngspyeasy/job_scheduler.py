#!/usr/bin/env python
import multiprocessing
import subprocess
import time
from threading import Thread
from Queue import Queue

from logger import logger
from job_dependency_tree import JobDependencyTree
import os

job_requests = Queue(-1)  # infinite shared job queue


class JobScheduler(Thread):
    def __init__(self, test_mode=False, timeout=60):
        super(JobScheduler, self).__init__()
        self.logger = logger()
        self.logger.debug("job_scheduler_init")

        numcores = multiprocessing.cpu_count()  # min=8
        self.logger.debug("numcores=%d", numcores)

        if numcores < 2:
            raise RuntimeError("Number of available cores %d (< 2).", numcores)

        numjobsallowed = numcores / 2
        self.logger.debug("[scheduler]: numjobsallowed=%d", numjobsallowed)

        self.processes = []
        self.tree = JobDependencyTree()
        self.max_processes = numjobsallowed
        self.test_mode = test_mode
        self.timeout = timeout
        self.all_done = False
        self.stop_all = False

    def run(self):
        global job_requests
        while not self.stop_all:
            while not job_requests.empty():
                (req_id, req_dep, req_command) = job_requests.get()
                job_requests.task_done()
                self.logger.debug("job_request_found: %s", req_id)

                if req_id == "stop_all":
                    self.logger.info("(stop_all) message received. Preparing to stop...")
                    self.stop_all = True
                    break

                if req_id == "all_done":
                    self.logger.info("(all_done) message received")
                    self.all_done = True

                try:
                    self.tree.append(req_id, req_dep, req_command)
                except ValueError, e:
                    self.logger.exception(e)
                    self.logger.info("[scheduler]: Preparing to stop...")
                    self.stop_all = True
                    break

            if len(self.processes) < self.max_processes:
                (job_id, job_command) = self.tree.get()

                if job_id is None and job_command is None:
                    if self.all_done and not self.tree.has_running_jobs():
                        self.logger.info("The job queue is empty. Preparing to stop..")
                        break

                elif job_command is None:
                    self.tree.done(job_id, 0)
                else:
                    self.logger.debug("job_to_run: %s", job_id)
                    self.logger.debug("command_to_run: [[\n%s \n]]", job_command)

                    # WARNING! using '-i' option with bash will create interactive shell which will have parent's stdin
                    # and get the SIGINT first..
                    if self.test_mode:
                        job_command = "pwd"
                    proc = subprocess.Popen(["/bin/bash", "-c", job_command], env=os.environ.copy())
                    self.processes.append((proc, job_command, job_id))

            self.update_processes()

        self.logger.info("stopping all running processes..")

        waiting_time = 0
        while len(self.processes) > 0:
            time.sleep(0.5)
            if waiting_time >= 2 * self.timeout:
                self.terminate_processes(sigkill=True)
            elif waiting_time >= self.timeout:
                self.terminate_processes()
            else:
                self.update_processes()
            waiting_time += 0.5

        self.logger.info("all stopped")

    def update_processes(self):
        unfinished = []
        for (p, c, job_id) in self.processes:
            ret = p.poll()
            if ret is None:
                unfinished.append((p, c, job_id))
            else:
                self.logger.debug("job_done: %s", job_id)
                self.tree.done(job_id, ret)
                if ret != 0:
                    self.logger.error("Command [[\n%s \n]] completed with error. See logs for details", c)

        self.processes = unfinished

    def terminate_processes(self, sigkill=False):
        unfinished = []
        for (p, c, job_id) in self.processes:
            ret = p.poll()
            if ret is None:
                unfinished.append((p, c, job_id))
                self.logger.error("Terminating process after timeout (%d) [[\n%s \n]]", self.timeout, c)
                if sigkill:
                    p.kill()
                else:
                    p.terminate()
        self.processes = unfinished


def submit(id, command=None, dependencies=None):
    global job_requests
    job_requests.put((id, dependencies, command))


def stop():
    submit("stop_all")


def all_done():
    submit("all_done")
