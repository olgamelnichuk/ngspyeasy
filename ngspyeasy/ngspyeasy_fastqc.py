#!/usr/bin/env python
import job_scheduler
from docker import docker_cmd
from settings import NGSEASYVERSION
from logger import log_info, log_set_current_step
import job_id_generator
import projects_dir


def ngspyeasy_fastqc(tsv_conf, projects_home, dependencies):
    log_set_current_step("ngspyeasy_fastqc")
    log_info("Schedule FastQC jobs")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        cmd = ["python /ngspyeasy/bin/ngspyeasy_fastqc_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        job_id = job_id_generator.get_next(["fastqc", sample_id])
        prev_job_ids = [x for x in [dependencies.get(sample_id, None)] if x is not None]

        log_info("New FastQC job(sample_id='%s', job_id='%s', dependencies='%s')" % (sample_id, job_id, prev_job_ids))

        job_scheduler.submit(
            job_id, docker_cmd(job_id, "compbio/ngseasy-fastqc:" + NGSEASYVERSION, " ".join(cmd), projects_home,
                               projects_dir.resources_dir(projects_home), pipeman=False),
            prev_job_ids)
        dependencies[sample_id] = job_id