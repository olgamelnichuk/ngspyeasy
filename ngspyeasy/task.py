class Task(object):
    def __init__(self):
        pass

    def execute(self, executor):
        samples2run = parallel_samples(task, vars)
        files2run = parallel_files(task, vars)

        if len(samples2run) > 0:
            for sample in samples2run:
                executor.submit(task_index, dict(curr_sample=sample))
        elif len(files2run) > 0:
            for file in files2run:
                executor.submit(task_index, dict(curr_file=file))
        else:
            executor.submit(task_index)

        executor.wait_for_jobs()
