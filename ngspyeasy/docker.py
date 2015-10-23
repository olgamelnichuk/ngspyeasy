#!/usr/bin/env python

USER = "pipeman"

HOME = "/home/" + USER

NGS_PROJECTS = HOME + "/ngs_projects"

NGS_RESOURCES = NGS_PROJECTS + "/ngseasy_resources"

DOCKER_OPTS = "-v /opt/ngspyeasy:/ngspyeasy:ro"


def docker_options(name, projects_home, resources_home, pipeman=True):
    options = ["--rm", "-P", "-w", HOME, "-e", "HOME=" + HOME]
    if pipeman:
        options.extend(["-e", "USER=" + USER, "--user", USER])

    options.extend(["--name", name])
    options.extend(["-v", projects_home + ":" + NGS_PROJECTS])
    options.extend(["-v", resources_home + ":" + NGS_RESOURCES])
    options.append(DOCKER_OPTS)
    return options


def wrap(name, image, cmd, projects_home, resources_home, pipeman=True):
    docker_run = ["docker", "run"] + docker_options(name, projects_home, resources_home, pipeman)
    docker_run.append(image)
    docker_run.append(cmd)
    return " ".join(docker_run)


def wrap_lsf(name, image, cmd, projects_home, resources_home, dependencies, pipeman=True):
    docker_image = "LSB_DOCKER_IMAGE=%s" % image
    docker_opts = "LSB_DOCKER_OPTIONS=\"%s\"" % " ".join(docker_options(name, projects_home, resources_home, pipeman))
    bsub_cmd = "bsub -J %s -w % -o out.log -e error.log %s" % (name, cmd, dependencies)
    return ";".join([docker_image, docker_opts, bsub_cmd])


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
        cmd = ["python /ngspyeasy/bin/%s" % self.executable, "-v" if self.verbose else "", "-c", self.config_name, "-d",
               NGS_PROJECTS, "-i", self.sample_id]

        if self.task:
            cmd += ["-t", self.task]
        return " ".join(cmd)