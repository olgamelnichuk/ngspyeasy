#!/usr/bin/env python

import job_scheduler
import projects_dir
from settings import NGSEASYVERSION
from docker import docker_cmd
from logger import log_set_current_step, log_info
from ngspyeasy import job_id_generator


def ngspyeasy_alignment(tsv_conf, projects_home, dependencies):
    log_set_current_step("ngspyeasy_alignment")
    log_info("Schedule alignment jobs")

    for row in tsv_conf.all_rows():
        sample_id = row.sample_id()
        aligner_type = row.aligner()

        cmd = ["python /ngspyeasy/bin/ngspyeasy_alignment_job.py", "-v", "-c", tsv_conf.filename(), "-d",
               "/home/pipeman/ngs_projects", "-i", sample_id]

        if aligner_type == "no-align":
            log_info("[%s] No alignment jobs to be run for sample: '%s'." % (aligner_type, sample_id))
            continue
        elif aligner_type == "bwa":
            commands = [("compbio/ngseasy-bwa:" + NGSEASYVERSION, cmd)]
        elif aligner_type == "novoalign":
            commands = [("compbio/ngseasy-novoalign:" + NGSEASYVERSION, cmd)]
        elif aligner_type == "stampy":
            commands = [
                ("compbio/ngseasy-bwa:" + NGSEASYVERSION, cmd + ["-p", "stampy_bwa"]),
                ("compbio/ngseasy-stampy:" + NGSEASYVERSION, cmd + ["-p", "stampy_stampy"]),
                ("compbio/ngseasy-picardtools:" + NGSEASYVERSION, cmd + ["-p", "stampy_picard1"]),
                ("compbio/ngseasy-picardtools:" + NGSEASYVERSION, cmd + ["-p", "stampy_picard2"])
            ]
        elif aligner_type == "bowtie2":
            commands = [("compbio/ngseasy-bowtie2:" + NGSEASYVERSION, cmd)]
        elif aligner_type == "snap":
            commands = [("compbio/ngseasy-snap:" + NGSEASYVERSION, cmd)]
        else:
            raise ValueError("Unknown aligner type: %s" % aligner_type)

        for (image, command) in commands:
            schedule(sample_id, projects_home, image, command, dependencies)


def schedule(sample_id, projects_home, docker_image, cmd, dependencies):
    job_id = job_id_generator.get_next(["aligner", sample_id])
    prev_job_ids = [x for x in [dependencies.get(sample_id, None)] if x is not None]

    log_info(
        "New aligner job(sample_id='%s', job_id='%s', dependencies=['%s'])" % (sample_id, job_id, prev_job_ids))

    job_scheduler.submit(
        job_id, docker_cmd(job_id, docker_image, " ".join(cmd), projects_home,
                           projects_dir.resources_dir(projects_home), pipeman=False), prev_job_ids)
    dependencies[sample_id] = job_id
