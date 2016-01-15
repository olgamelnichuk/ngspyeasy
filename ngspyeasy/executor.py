import multiprocessing
import subprocess
import time
import sys
import traceback

import re
from logger import logger
import os


class TaskSubmitError(Exception):
    pass


class Provider(object):
    def submit(self, name, cmd):
        raise NotImplementedError()

    def list(self):
        raise NotImplementedError()

    def stop(self):
        raise NotImplementedError()


class LSFProvider(Provider):
    def __init__(self, queue=None, log_dir=None):
        self._queue = queue
        self._log_dir = log_dir

    def submit(self, name, cmd):
        command = ["bsub"]

        command.extend(["-J", name])

        if self._queue:
            command.extend(["-q", str(self._queue)])

        if self._log_dir:
            stderr = os.path.join(self._log_dir, "lsf-%J.err")
            stdout = os.path.join(self._log_dir, "lsf-%J.out")
            command.extend(["-o", stdout])
            command.extend(["-e", stderr])

        command.append(cmd)
        logger().debug("Submitting job with :%s %s", command)

        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = "".join([l for l in proc.stdout])
        err = "".join([l for l in proc.stderr])
        expr = 'Job <(?P<job_id>.+)> is submitted.*'
        match = re.search(expr, out)
        if proc.wait() != 0 or not match:
            raise TaskSubmitError("%s\n"
                                  "Executed command:\n%s\n%s\n" % (
                                      out,
                                      err,
                                      " ".join(cmd)
                                  ))
        return match.group('job_id')

    def list(self):
        proc = subprocess.Popen([], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        jobs = []
        for l in proc.stdout:
            fields = [x for x in l.strip().split(" ") if x]
            try:
                long(fields[0])
            except:
                continue
            jobs.append(fields[0])
        err = "".join([l for l in proc.stderr])
        if proc.wait() != 0:
            raise ValueError("Error while listing jobs:\n%s" % (err))
        return jobs

    def stop(self):
        jobs = self.list()
        for job_id in jobs:
            self._cancel(job_id)

    def _cancel(self, job_id):
        cmd = ["bkill", str(job_id)]
        subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


class LocalProvider(Provider):
    def __init__(self):
        self._pool_size = max(multiprocessing.cpu_count() / 2, 1)
        self._procs = []
        self._waiting_jobs = []

    def submit(self, name, cmd):
        self._waiting_jobs.insert(0, (name, cmd))
        self._run_next()
        return name

    def list(self):
        self._run_next()
        return [x[0] for x in self._waiting_jobs] + [x[2] for x in self._procs]

    def stop(self):
        while (len(self._procs) > 0):
            self._stop_all()
            time.sleep(0.5)

    def _stop_all(self, sigkill=False):
        self._update()
        for (proc, cmd, name) in self._procs:
            if sigkill:
                proc.kill()
            else:
                proc.terminate()

    def _run_next(self):
        self._update()
        if len(self._procs) > self._pool_size:
            return
        if len(self._waiting_jobs) == 0:
            return
        (name, cmd) = self._waiting_jobs.pop()
        proc = subprocess.Popen(["/bin/bash", "-c", cmd], env=os.environ.copy())
        self._procs.append((proc, cmd, name))

    def _update(self):
        unfinished = [x for x in self._procs if x[0].poll() is None]
        self._procs = unfinished


work_queue = multiprocessing.Queue(-1)

results_queue = multiprocessing.Queue(-1)


def start(provider, log_dir):
    e = JobExecutor(provider=provider, log_dir=log_dir)
    e.start()


def stop():
    submit("STOP", None)


def submit(name, cmd):
    work_queue.put((name, cmd))


class JobExecutor(multiprocessing.Process):
    def __init__(self, provider, log_dir):
        super(JobExecutor, self).__init__()
        self._provider = LSFProvider(log_dir=log_dir) if provider == "lsf" else LocalProvider()
        self._running = True
        self._mapping = dict()
        self._running_jobs = []

    def run(self):
        try:
            self.run_with_exception()
        except:
            results_queue.put("STOP")
            (type, value, tb) = sys.exc_info()
            e = "".join(traceback.format_exception(type, value, tb))
            logger().error("executor: exiting with exception:\n %s" % e)

    def run_with_exception(self):
        global work_queue
        global results_queue

        while self._running:
            if not work_queue.empty():
                (name, cmd) = work_queue.get(block=False)

                if name == "STOP":
                    logger().info("executor: received [STOP] message")
                    results_queue.put("STOP")
                    self._provider.stop()
                    break

                logger().info("executor: received cmd to run: name=%s" % name)
                self._submit(name, cmd)
            self._update_results()

    def _submit(self, name, cmd):
        job_id = self._provider.submit(name, cmd)
        self._mapping[job_id] = name
        self._running_jobs.append(job_id)

    def _update_results(self):
        job_list = self._provider.list()
        finished_jobs = [x for x in self._running_jobs not in set(job_list)]
        for job_id in finished_jobs:
            results_queue.put(self._mapping[job_id])
        self._running_jobs = job_list
