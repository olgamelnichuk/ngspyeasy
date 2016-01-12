import job_id_generator
import jinja2
import glob


# TODO

def as_commands(task_index, task, vars):
    samples2run = parallel_samples(task, vars)
    files2run = parallel_files(task, vars)


def parallel_samples(task, vars):
    tmpl = task.get("samples", None)
    if tmpl is None:
        return []
    samples_str = jinja2.Template(tmpl).render(vars)
    return eval(samples_str)


def parallel_files(task, vars):
    tmpl = task.get("files", None)
    if tmpl is None:
        return []
    pattern = jinja2.Template(tmpl).render(vars)
    return sorted(glob.glob(pattern))


def job_id(task_index):
    return job_id_generator.get_next(["task_" + str(task_index)])


def cmd(task_index, run_index, pipeline_script, options):
    executable = "ngspyeasy_task"
    return " ".join([executable,
                     pipeline_script,
                     "--task_index", str(task_index),
                     "--run_index"
                     ]) + options
