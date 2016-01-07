#!/usr/bin/env python

###
# Copyright 2015, EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###

import shutil

import shcmd
import os
from utils import uniq_set
from settings import VERSION
from logger import logger


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

    def log_path(self, name):
        return os.path.join(self.log_dir(), name)

    def log_dir(self):
        return os.path.join(self.projects_home, "run_log")

    def config_dir(self):
        return os.path.join(self.projects_home, "config_files")

    def raw_fastq_dir(self):
        return os.path.join(self.projects_home, "raw_fastq")

    def tmp_dir(self):
        return os.path.join(self.projects_home, "tmp")

    def raw_fastq_path(self, fastq_name):
        return os.path.join(self.raw_fastq_dir(), fastq_name)

    def resources_dir(self):
        return os.path.join(self.projects_home,
                            "ngseasy_resources") if self.resources_home is None else self.resources_home

    def project_dir(self, project_id):
        return os.path.join(self.projects_home, project_id)

    def sample_dir(self, project_id, sample_id):
        return os.path.join(self.project_dir(project_id), sample_id)

    def fix_file_permissions(self, project_id, sample_id, uid, gid):
        dir = self.sample_dir(project_id, sample_id)
        logger().info("chown -r %s %s %s" % (str(uid), str(gid), dir))
        shcmd.chown(dir, uid, gid)

    def init_structure(self, tsv_conf):
        uniq_sample_dirs = uniq_set(
            [self.sample_dir(x.project_id(), x.sample_id()) for x in tsv_conf.all_rows()])

        logger().info("Checking the project directory structure...")
        logger().debug("Sample dirs to be checked:\n %s" % "\n".join(uniq_sample_dirs))

        def makedir_ifnotexist(dir):
            if not os.path.isdir(dir):
                logger().info("Creating directory: %s" % dir)
                os.makedirs(dir)

        makedir_ifnotexist(self.tmp_dir())

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

    def check_fastq(self, tsv_conf):
        logger().info("Checking if we need to move raw FastQ files...")
        for row in tsv_conf.all_rows():
            sample_dir = SampleDir(self.sample_dir(row.project_id(), row.sample_id()))

            for fq in [row.fastq1(), row.fastq2()]:
                move_fastq(self.raw_fastq_path(fq), sample_dir.fastq_path(fq))


def move_fastq(source, dest):
    if not os.path.isfile(dest):
        logger().info("FastQ does not exist: %s. Checking if the raw FastQ exist..." % dest)
        if not os.path.isfile(source):
            raise IOError("FastQ file does not exist: %s" % source)

        logger().info("Moving FastQ file: %s --> %s" % (source, dest))
        shutil.move(source, dest)

    logger().info("OK (FastQ file exists: %s)" % dest)


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
