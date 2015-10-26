#!/usr/bin/env python
import os

USER = "pipeman"

HOME = "/home/" + USER

NGS_PROJECTS = HOME + "/ngs_projects"

NGS_RESOURCES = NGS_PROJECTS + "/ngseasy_resources"

DOCKER_OPTS = ""


def wrap(name, image, cmd, projects_home, resources_home, pipeman=True):
    docker_run = ["docker", "run", "--rm", "-P", "-w", HOME, "-e", "HOME=" + HOME]
    if pipeman:
        docker_run.extend(["-e", "USER=" + USER, "--user", USER])

    ngspyeasy_home = os.path.dirname(os.path.realpath(__file__))
    docker_run.extend(["--name", name])
    docker_run.extend(["-v", projects_home + ":" + NGS_PROJECTS])
    docker_run.extend(["-v", resources_home + ":" + NGS_RESOURCES])
    docker_run.extend(["-v", ngspyeasy_home + ":/ngspyeasy:ro"])
    docker_run.append(DOCKER_OPTS)
    docker_run.append(image)
    docker_run.append(cmd)
    return " ".join(docker_run)


class JobCommand(object):
    def __init__(self, executable, config_name, sample_id, **kwargs):
        self.executable = executable
        self.sample_id = sample_id
        self.config_name = config_name
        self.verbose = kwargs.get("verbose", False)
        self.task = kwargs.get("task", None)

    def add_task(self, task):
        return JobCommand(self.executable, self.config_name, self.sample_id, verbose=self.verbose, task=task)

    def as_string(self):
        cmd = ["python /ngspyeasy/%s" % self.executable, "-v" if self.verbose else "", "-c", self.config_name, "-d",
               NGS_PROJECTS, "-i", self.sample_id]

        if self.task:
            cmd += ["-t", self.task]
        return " ".join(cmd)