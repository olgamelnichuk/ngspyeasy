#!/usr/bin/env python

import os

from logger import log_debug
from utils import uniq_set


def get_log_dir(projects_home):
    return os.path.join(projects_home, "run_log")


def get_config_dir(projects_home):
    return os.path.join(projects_home, "config_files")


def get_config_path(projects_home, tsv_name):
    return os.path.join(get_config_dir(projects_home), tsv_name)


def get_project_dir(projects_home, project_id):
    return os.path.join(projects_home, project_id)


def get_sample_dir(projects_home, project_id, sample_id):
    return os.path.join(get_project_dir(projects_home, project_id), sample_id)


def get_sample_fastq_dir(sample_dir):
    return os.path.join(sample_dir, "fastq")


def get_sample_tmp_dir(sample_dir):
    return os.path.join(sample_dir, "tmp")


def get_sample_alignments_dir(sample_dir):
    return os.path.join(sample_dir, "alignments")


def get_sample_vcf_dir(sample_dir):
    return os.path.join(sample_dir, "vcf")


def get_sample_reports_dir(sample_dir):
    return os.path.join(sample_dir, "reports")


def get_sample_config_dir(sample_dir):
    return os.path.join(sample_dir, "config_files")


def get_sample_log_dir(sample_dir):
    return os.path.join(sample_dir, "run_logs")


def makedir_ifnotexist(dir):
    if not os.path.isdir(dir):
        log_debug("Creating dir: %s", dir)
        os.makedirs(dir)


def init_project(tsv_conf, projects_home):
    uniq_sample_dirs = uniq_set(
        map(lambda x: get_sample_dir(projects_home, x.get_project_id(), x.get_sample_id()), tsv_conf.all_rows()))

    log_debug("unique sample dirs to check/create: %s", "\n".join(uniq_sample_dirs))

    for sample_dir in uniq_sample_dirs:
        makedir_ifnotexist(sample_dir)
        makedir_ifnotexist(get_sample_fastq_dir(sample_dir))
        makedir_ifnotexist(get_sample_tmp_dir(sample_dir))
        makedir_ifnotexist(get_sample_alignments_dir(sample_dir))
        makedir_ifnotexist(get_sample_vcf_dir(sample_dir))
        makedir_ifnotexist(get_sample_reports_dir(sample_dir))
        makedir_ifnotexist(get_sample_config_dir(sample_dir))
        makedir_ifnotexist(get_sample_log_dir(sample_dir))

    for r,d,f in os.walk(projects_home):
        os.chmod(r, 0775)