#!/usr/bin/env python
import subprocess
import sys

from logger import log_debug, log_error, log_info
import os.path
from projects_dir import config_full_path, config_dir


def check_cmdline_options(tsv_config_file, ngs_projects_dir):
    (projects_home, errmsg) = check_ngs_projects_dir_option(ngs_projects_dir)
    if errmsg is not None:
        return None, None, errmsg

    (tsv_name, errmsg) = check_tsv_config_file_option(tsv_config_file, projects_home)
    if errmsg is not None:
        return None, None, errmsg

    return tsv_name, projects_home, None


def check_ngs_projects_dir_option(ngs_projects_dir):
    if not ngs_projects_dir:
        return None, "Projects directory is a required parameter. Can't find one."

    projects_home = os.path.abspath(ngs_projects_dir)

    if not os.path.isdir(projects_home):
        return None, "Projects directory '" + projects_home + "' does not exist."

    return projects_home, None


def check_tsv_config_file_option(tsv_config_file, projects_home):
    if not tsv_config_file:
        return None, None, "TSV config file is not specified."

    expected_path = config_full_path(projects_home, os.path.basename(tsv_config_file))

    if os.path.isabs(tsv_config_file):
        if os.path.abspath(tsv_config_file) != expected_path:
            return None, "Config file must be in the projects config directory: '" + config_dir(
                projects_home) + "'"

    if not os.path.isfile(expected_path):
        return None, "Config file '" + expected_path + "' does not exist."

    return os.path.basename(expected_path), None


def run_command(cmd):
    proc = subprocess.Popen(["/bin/bash", "-i", "-c", "source ~/.bashrc; " + " ".join(cmd)],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)

    stdout = []
    try:
        for line in iter(proc.stdout.readline, ''):
            sys.stdout.write(line)
            sys.stdout.flush()
            stdout.append(line)
    except KeyboardInterrupt:
        log_info("got Ctrl + C")

    log_debug("cmd: \n" + "".join(stdout))

    if proc.returncode:
        log_error("Command [[\n%s\n]] failed. See logs for details", " ".join(cmd))
