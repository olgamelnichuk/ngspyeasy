#!/usr/bin/env python

import shutil

import shcmd
import os
from utils import uniq_set
from settings import VERSION
from logger import log_info, log_debug


class ProjectsDir(object):
    def __init__(self, projects_home, resources_home=None):
        self.projects_home = projects_home
        self.resources_home = resources_home

    def root(self):
        return self.projects_home

    def sample_log_file(self, tsv_name, sample_id="ALL_SAMPLES"):
        log_name = tsv_name + '@' + str(sample_id)
        return os.path.join(self.log_dir(), 'ngspyeasy.' + VERSION + '.' + log_name + '.log')

    def main_log_file(self, tsv_name):
        return os.path.join(self.log_dir(), 'ngspyeasy.' + VERSION + '.' + tsv_name + '.log')

    def config_path(self, tsv_name):
        return os.path.join(self.config_dir(), tsv_name)

    def resource_path(self, name):
        return os.path.join(self.resources_dir(), name)

    def log_dir(self):
        return os.path.join(self.projects_home, "run_log")

    def config_dir(self):
        return os.path.join(self.projects_home, "config_files")

    def raw_fastq_dir(self):
        return os.path.join(self.projects_home, "raw_fastq")

    def raw_fastq_path(self, fastq_name):
        return os.path.join(self.raw_fastq_dir(), fastq_name)

    def resources_dir(self):
        return os.path.join(self.projects_home,
                            "ngseasy_resources") if self.resources_home is None else self.resources_home

    def project_dir(self, project_id):
        return os.path.join(self.projects_home, project_id)

    def sample_dir(self, project_id, sample_id):
        return os.path.join(self.project_dir(project_id), sample_id)

    def fix_file_permissions(self, project_id, sample_id, logger):
        dir = self.sample_dir(project_id, sample_id)
        logger.info("chmod -r 0777 0666 %s" % dir)
        shcmd.chmod(dir, 0777, 0666)

    def init_structure(self, tsv_conf):
        uniq_sample_dirs = uniq_set(
            [self.sample_dir(x.project_id(), x.sample_id()) for x in tsv_conf.all_rows()])

        log_info("Checking the project directory structure...")
        log_debug("Sample dirs to be checked:\n %s" % "\n".join(uniq_sample_dirs))

        def makedir_ifnotexist(dir):
            if not os.path.isdir(dir):
                log_info("Creating directory: %s" % dir)
                os.makedirs(dir)

        for dir in uniq_sample_dirs:
            sample_dir = SampleDir(dir)
            makedir_ifnotexist(sample_dir.root())
            makedir_ifnotexist(sample_dir.fastq_dir())
            makedir_ifnotexist(sample_dir.tmp_dir())
            makedir_ifnotexist(sample_dir.alignments_dir())
            makedir_ifnotexist(sample_dir.vcf_dir())
            makedir_ifnotexist(sample_dir.reports_dir())
            makedir_ifnotexist(sample_dir.config_dir())
            makedir_ifnotexist(sample_dir.log_dir())

        project_dirs = [self.project_dir(x) for x in uniq_set([y.project_id() for y in tsv_conf.all_rows()])]

        for d in project_dirs:
            log_info("Chmod 0775 on everything under %s" % d)
            shcmd.chmod(d, 0775, 0664)

    def check_fastq(self, tsv_conf):
        log_info("Checking if we need to move raw FastQ files...")
        for row in tsv_conf.all_rows():
            sample_dir = SampleDir(self.sample_dir(row.project_id(), row.sample_id()))

            for fq in [row.fastq1(), row.fastq2()]:
                move_fastq(self.raw_fastq_path(fq), sample_dir.fastq_path(fq))


def move_fastq(source, dest):
    if not os.path.isfile(dest):
        log_info("FastQ does not exist: %s. Checking if the raw FastQ exist..." % dest)
        if not os.path.isfile(source):
            raise IOError("FastQ file does not exist: %s" % source)

        log_info("Moving FastQ file: %s --> %s" % (source, dest))
        shutil.move(source, dest)

    log_info("OK (FastQ file exists: %s)" % dest)


class SampleDir(object):
    def __init__(self, sample_home):
        self.sample_home = sample_home

    def root(self):
        return self.sample_home

    def fastq_dir(self):
        return os.path.join(self.sample_home, "fastq")

    def tmp_dir(self):
        return os.path.join(self.sample_home, "tmp")

    def alignments_dir(self):
        return os.path.join(self.sample_home, "alignments")

    def vcf_dir(self):
        return os.path.join(self.sample_home, "vcf")

    def reports_dir(self):
        return os.path.join(self.sample_home, "reports")

    def config_dir(self):
        return os.path.join(self.sample_home, "config_files")

    def log_dir(self):
        return os.path.join(self.sample_home, "run_logs")

    def fastq_path(self, filename):
        return os.path.join(self.fastq_dir(), filename)

    def alignments_path(self, filename):
        return os.path.join(self.alignments_dir(), filename)

    def reports_path(self, filename):
        return os.path.join(self.reports_dir(), filename)

    def vcf_path(self, filename):
        return os.path.join(self.vcf_dir(), filename)
