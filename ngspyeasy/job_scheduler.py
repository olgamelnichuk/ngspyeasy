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

job_requests = Queue(-1)  # infinite shared job queue


class JobScheduler(Thread):
    def __init__(self, logfile=None):
        super(JobScheduler, self).__init__()
        self.logger = self.create_logger(logfile)
        self.logger.debug("job_scheduler_init")

        numcores = multiprocessing.cpu_count()  # min=8
        self.logger.debug("numcores=%d", numcores)

        if numcores < 2:
            raise RuntimeError("Number of available cores %d (< 2).", numcores)

        numjobsallowed = numcores / 2
        self.logger.debug("numjobsallowed=%d", numjobsallowed)

        self.processes = []
        self.tree = JobDependencyTree()
        self.max_processes = numjobsallowed
        self.stopped = False

    def create_logger(self, logfile):
        logger = logging.getLogger("JobScheduler")
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(threadName)s %(levelname)s [JobScheduler] %(message)s')

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
                self.logger.debug("job_request_found: %s", req_id)

                if req_id == "stop_all":
                    pass
                    #self.logger.info("[stop_all] message received. Preparing to stop...")
                    #self.stopped = True
                    #break
                else:
                    try:
                        self.tree.append(req_id, req_dep, req_command)
                    except ValueError as ex:
                        self.logger.error(ex.message, ex.args)
                        self.logger.info("Preparing to stop...")
                        self.stopped = True
                        break

            while len(self.processes) >= self.max_processes:
                time.sleep(0.2)
                self.update_processes()

            (job_id, job_command) = self.tree.get()
            if job_id is None or job_command is None:
                continue

            self.logger.debug("job_to_run: %s", job_id)
            self.logger.debug("command_to_run: [[\n %s \n]]", job_command)

            proc = subprocess.Popen(shlex.split(job_command))
            self.processes.append((proc, job_command, job_id))

        self.logger.info("[stop_all] waiting processes to finish..")
        while len(self.processes) > 0:
            time.sleep(0.2)
            self.update_processes()

        self.logger.info("[stop_all] all stopped")

    def update_processes(self):
        unfinished = []
        for p, c, job_id in self.processes:
            ret = p.poll()
            if ret is None:
                unfinished.append((p, c))
            else:
                self.logger.debug("job_done: %s", job_id)
                self.tree.done(job_id, ret)
                if ret != 0:
                    self.logger.error("Command (%s) completed with error. See logs for details", c)

        self.processes = unfinished


def submit(id, command, dependencies=None):
    global job_requests
    log_debug("job_request_submit: (%s, %s)", id, command)
    job_requests.put((id, dependencies, command))
    log_debug("job_requests_size: %d", job_requests.qsize())


def stop():
    submit("stop_all", None, None)
