from threading import Thread
from Queue import Queue
from ngspyeasy.providers import LocalProvider


class TaskExecutor(Thread):
    def __init__(self, provider):
        super(TaskExecutor, self).__init__()
        self._provider = LocalProvider()
        self._requests = Queue()
        self._running = True
        self._id_map = dict()
        self._running_jobs = []

    def submit(self, name, cmd):
        self._requests.put((name, cmd))

    def terminate(self):
        self._requests.put(("terminate", None))

    def wait_for_jobs(self):
        self._requests.join()

    def run(self):
        while self._running:
            while self._requests.empty():
                (name, cmd) = self._requests.get()

                if name == "terminate":
                    self._running = False
                    break

                job_id = self._provider.submit(name, cmd)
                self._id_map[name] = job_id
                self._running_jobs.append(job_id)

            job_list = self._provider.list()
            diff = [x for x in self._running_jobs not in set(job_list)]
            for j in diff:
                self._requests.task_done()
            self._running_jobs = job_list

        self._provider.terminate()
