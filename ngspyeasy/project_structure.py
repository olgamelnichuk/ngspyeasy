#!/usr/bin/env python

import os
from shutil import move

from logger import log_debug, log_info
from utils import uniq_set


def get_log_dir(projects_home):
    return os.path.join(projects_home, "run_log")


def get_config_dir(projects_home):
    return os.path.join(projects_home, "config_files")


def get_raw_fastq_dir(projects_home):
    return os.path.join(projects_home, "raw_fastq")


def get_resources_dir(projects_home):
    return os.path.join(projects_home, "ngseasy_resources")


def get_config_path(projects_home, tsv_name):
    return os.path.join(get_config_dir(projects_home), tsv_name)


def get_project_dir(projects_home, project_id):
    return os.path.join(projects_home, project_id)


def get_sample_dir(projects_home, project_id, sample_id):
    return os.path.join(get_project_dir(projects_home, project_id), sample_id)


def get_sample_fastq_dir(sample_dir):
    return os.path.join(sample_dir, "fastq")


def get_sample_fastq_path(sample_dir, fastq_filename):
    return os.path.join(get_sample_fastq_dir(sample_dir), fastq_filename)


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

    log_debug("Unique sample dirs to check/create: %s", "\n".join(uniq_sample_dirs))
    log_info("Creating required project structure...")

    for sample_dir in uniq_sample_dirs:
        makedir_ifnotexist(sample_dir)
        makedir_ifnotexist(get_sample_fastq_dir(sample_dir))
        makedir_ifnotexist(get_sample_tmp_dir(sample_dir))
        makedir_ifnotexist(get_sample_alignments_dir(sample_dir))
        makedir_ifnotexist(get_sample_vcf_dir(sample_dir))
        makedir_ifnotexist(get_sample_reports_dir(sample_dir))
        makedir_ifnotexist(get_sample_config_dir(sample_dir))
        makedir_ifnotexist(get_sample_log_dir(sample_dir))

    log_info("Chmod 0775 on everything under %s", projects_home)

    for r, d, f in os.walk(projects_home):
        os.chmod(r, 0775)


def init_fastq(tsv_conf, projects_home):
    for row in tsv_conf.all_rows():
        fastq1 = row.get_fastq1()
        fastq2 = row.get_fastq2()

        sample_home = get_sample_dir(projects_home, row.get_project_id(), row.get_sample_id())

        move_fastq(projects_home, sample_home, fastq1)
        move_fastq(projects_home, sample_home, fastq2)


def move_fastq(projects_home, sample_home, fastq_name):
    source = os.path.join(get_raw_fastq_dir(projects_home), fastq_name)
    dest = os.path.join(get_sample_fastq_dir(sample_home), fastq_name)

    log_info("Checking if fastq file exists: %s", dest)

    if not os.path.isfile(dest):
        if not os.path.isfile(source):
            raise IOError("Fastq file doesn't exist: " + source)

        log_info("Moving fastq file: src=%s, dest=%s", source, dest)
        move(source, dest)
    else:
        log_info("OK (fastq file exists: %s)", dest)
