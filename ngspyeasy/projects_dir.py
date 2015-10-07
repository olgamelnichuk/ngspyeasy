#!/usr/bin/env python

import os
import shutil

from logger import log_debug, log_info
from utils import uniq_set
from settings import VERSION


def log_dir(projects_home):
    return os.path.join(projects_home, "run_log")


def config_dir(projects_home):
    return os.path.join(projects_home, "config_files")


def raw_fastq_dir(projects_home):
    return os.path.join(projects_home, "raw_fastq")


def resources_dir(projects_home):
    return os.path.join(projects_home, "ngseasy_resources")


def project_dir(projects_home, project_id):
    return os.path.join(projects_home, project_id)


def sample_dir(projects_home, project_id, sample_id):
    return os.path.join(project_dir(projects_home, project_id), sample_id)


def sample_fastq_dir(sample_dir):
    return os.path.join(sample_dir, "fastq")


def sample_tmp_dir(sample_dir):
    return os.path.join(sample_dir, "tmp")


def sample_alignments_dir(sample_dir):
    return os.path.join(sample_dir, "alignments")


def sample_vcf_dir(sample_dir):
    return os.path.join(sample_dir, "vcf")


def sample_reports_dir(sample_dir):
    return os.path.join(sample_dir, "reports")


def sample_config_dir(sample_dir):
    return os.path.join(sample_dir, "config_files")


def sample_log_dir(sample_dir):
    return os.path.join(sample_dir, "run_logs")


def sample_log_file(projects_home, tsv_name, sample_id):
    log_name = tsv_name + "@" + str(sample_id)
    return os.path.join(log_dir(projects_home), "ngspyeasy." + VERSION + "." + log_name + ".log")


def main_log_file(projects_home, tsv_name):
    return os.path.join(log_dir(projects_home), "ngspyeasy." + VERSION + "." + tsv_name + ".log")


def config_full_path(projects_home, tsv_name):
    return os.path.join(config_dir(projects_home), tsv_name)


def fastq_path(sample_dir, filename):
    return os.path.join(sample_fastq_dir(sample_dir), filename)


def makedir_ifnotexist(dir):
    if not os.path.isdir(dir):
        log_debug("Creating dir: %s", dir)
        os.makedirs(dir)


def init(projects_home, tsv_conf):
    uniq_sample_dirs = uniq_set(
        map(lambda x: sample_dir(projects_home, x.project_id(), x.sample_id()), tsv_conf.all_rows()))

    log_debug("Unique sample dirs to check/create: %s", "\n".join(uniq_sample_dirs))
    log_info("Creating required project structure...")

    for dir in uniq_sample_dirs:
        makedir_ifnotexist(dir)
        makedir_ifnotexist(sample_fastq_dir(dir))
        makedir_ifnotexist(sample_tmp_dir(dir))
        makedir_ifnotexist(sample_alignments_dir(dir))
        makedir_ifnotexist(sample_vcf_dir(dir))
        makedir_ifnotexist(sample_reports_dir(dir))
        makedir_ifnotexist(sample_config_dir(dir))
        makedir_ifnotexist(sample_log_dir(dir))

    log_info("Chmod 0775 on everything under %s", projects_home)

    for r, d, f in os.walk(projects_home):
        os.chmod(r, 0775)


def init_fastq(tsv_conf, projects_home):
    for row in tsv_conf.all_rows():
        fastq1 = row.fastq1()
        fastq2 = row.fastq2()

        sample_home = sample_dir(projects_home, row.project_id(), row.sample_id())

        move_fastq(projects_home, sample_home, fastq1)
        move_fastq(projects_home, sample_home, fastq2)


def move_fastq(projects_home, sample_home, fastq_name):
    source = os.path.join(raw_fastq_dir(projects_home), fastq_name)
    dest = os.path.join(sample_fastq_dir(sample_home), fastq_name)

    log_info("Checking if fastq file exists: %s", dest)

    if not os.path.isfile(dest):
        if not os.path.isfile(source):
            raise IOError("Fastq file doesn't exist: " + source)

        log_info("Moving fastq file: %s --> %s", source, dest)
        shutil.move(source, dest)
    else:
        log_info("OK (fastq file exists: %s)", dest)
