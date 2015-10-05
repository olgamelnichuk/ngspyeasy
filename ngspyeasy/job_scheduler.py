#!/usr/bin/env python
import sys
import multiprocessing
import shlex
import subprocess
import time
import logging
from threading import Thread

from Queue import Queue
from logger import log_debug
from job_dependency_tree import JobDependencyTree
import os

job_requests = Queue(-1)  # infinite shared job queue


class JobScheduler(Thread):
    def __init__(self, logfile=None, timeout=60):
        super(JobScheduler, self).__init__()
        self.logger = self.create_logger(logfile)
        self.logger.debug("[scheduler]: job_scheduler_init")

        numcores = multiprocessing.cpu_count()  # min=8
        self.logger.debug("[scheduler]: numcores=%d", numcores)

        if numcores < 2:
            raise RuntimeError("Number of available cores %d (< 2).", numcores)

        numjobsallowed = numcores / 2
        self.logger.debug("[scheduler]: numjobsallowed=%d", numjobsallowed)

        self.processes = []
        self.tree = JobDependencyTree()
        self.max_processes = numjobsallowed
        self.stopped = False
        self.timeout = timeout

    def create_logger(self, logfile):
        logger = logging.getLogger("JobScheduler")
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(threadName)s %(levelname)s %(message)s')

        if logfile:
            file_handler = logging.FileHandler(logfile, mode='w')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        else:
            stream_handler = logging.StreamHandler(sys.stdout)
            stream_handler.setFormatter(formatter)
            logger.addHandler(stream_handler)

        return logger

    def run(self):
        global job_requests
        while not self.stopped:
            while not job_requests.empty():
                (req_id, req_dep, req_command) = job_requests.get()
                job_requests.task_done()
                self.logger.debug("[scheduler]: job_request_found: %s", req_id)

                if req_id == "stop_all":
                    self.logger.info("[scheduler]: (stop_all) message received. Preparing to stop...")
                    self.stopped = True
                    break

                try:
                    self.tree.append(req_id, req_dep, req_command)
                except ValueError, e:
                    self.logger.exception(e)
                    self.logger.info("[scheduler]: Preparing to stop...")
                    self.stopped = True
                    break

            if len(self.processes) < self.max_processes:
                (job_id, job_command) = self.tree.get()
                if job_command is None:
                    if job_id is not None:
                        self.tree.done(job_id, 0)
                else:
                    self.logger.debug("[scheduler]: job_to_run: %s", job_id)
                    self.logger.debug("[scheduler]: command_to_run: [[\n%s \n]]", job_command)

                    proc = subprocess.Popen(["/bin/bash", "-i", "-c", job_command], env=os.environ.copy())
                    self.processes.append((proc, job_command, job_id))

            self.update_processes()

        self.logger.info("[scheduler]: stopping all running processes..")

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

        self.logger.info("[scheduler]: all stopped")

    def update_processes(self):
        unfinished = []
        for (p, c, job_id) in self.processes:
            ret = p.poll()
            if ret is None:
                unfinished.append((p, c, job_id))
            else:
                self.logger.debug("[scheduler]: job_done: %s", job_id)
                self.tree.done(job_id, ret)
                if ret != 0:
                    self.logger.error("[scheduler]: Command [[\n%s \n]] completed with error. See logs for details", c)

        self.processes = unfinished

    def terminate_processes(self, sigkill=False):
        unfinished = []
        for (p, c, job_id) in self.processes:
            ret = p.poll()
            if ret is None:
                unfinished.append((p, c, job_id))
                self.logger.error("[scheduler]: Terminating process after timeout (%d) [[\n%s \n]]", self.timeout, c)
                if sigkill:
                    p.kill()
                else:
                    p.terminate()
        self.processes = unfinished


def submit(id, command=None, dependencies=None):
    global job_requests
    log_debug("job_request_submit: (%s, %s)", id, command)
    job_requests.put((id, dependencies, command))
    log_debug("job_requests_size: %d", job_requests.qsize())


def stop():
    submit("stop_all")
