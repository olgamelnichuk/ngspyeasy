#!/usr/bin/env python

import projects_dir
import job_id_generator
import job_scheduler
from docker import docker_cmd
from logger import log_set_current_step, log_info
from settings import NGSEASYVERSION


def ngspyeasy_trimmomatic(tsv_conf, projects_home, dependencies):
    log_set_current_step("ngspyeasy_trimmomatic")
    log_info("Start: Trimmomatic")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        trim_type = row.trim()

        if trim_type not in ["atrim", "btrim"]:
            raise ValueError("Unknown trimmomatic type: %s" % trim_type)

        # trimmomatic step 1

        cmd = ["python /ngspyeasy/bin/ngspyeasy_trimmomatic_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        job_id = job_id_generator.get_next(["trimmomatic", sample_id])
        prev_job_ids = [dependencies[sample_id]]

        job_scheduler.submit(
            job_id, docker_cmd(job_id, "compbio/ngseasy-trimmomatic:" + NGSEASYVERSION, " ".join(cmd), projects_home,
                               projects_dir.resources_dir(projects_home), pipeman=False),
            [x for x in prev_job_ids if x is not None])
        dependencies[sample_id] = job_id

        # trimmomatic step 2

        cmd = ["python /ngspyeasy/bin/ngspyeasy_trimmomatic_fastqc_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        job_id = job_id_generator.get_next(["trimmomatic_fastqc", sample_id])
        prev_job_ids = [dependencies[sample_id]]

        job_scheduler.submit(
            job_id, docker_cmd(job_id, "compbio/ngseasy-fastqc:" + NGSEASYVERSION, " ".join(cmd), projects_home,
                               projects_dir.resources_dir(projects_home), pipeman=False),
            [x for x in prev_job_ids if x is not None])
        dependencies[sample_id] = job_id
