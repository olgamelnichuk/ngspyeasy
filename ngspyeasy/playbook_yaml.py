import copy
import glob

import os

import job_id_generator
import jinja2
from logger import logger
import yaml
import tsv_config


class PlayBookYaml(object):
    def __init__(self, cmd, plays, samples, variables):
        self._cmd = cmd
        self._plays = plays
        self._samples = samples
        self._vars = variables

    def plays(self):
        for index, play in enumerate(self._plays, start=0):
            yield self._create_play(index, play)

    def play_run(self, play_index, run_index):
        yaml_obj = self._plays[play_index]
        play = self._create_play(play_index, yaml_obj)
        play_run = play.play_run(run_index)
        return play_run.vars(self._vars), play_run.yaml(yaml_obj), play_run.name()

    def _create_play(self, index, yaml_obj):
        files = self._files2run(yaml_obj, self._vars)
        samples = self._samples2run(yaml_obj, self._vars)
        return PlayYaml(index, files, samples, self._cmd)

    @staticmethod
    def _samples2run(yaml_obj, variables):
        tmpl = yaml_obj.get("samples", None)
        if tmpl is None:
            return []
        samples_str = jinja2.Template(tmpl).render(variables)
        return eval(samples_str)

    @staticmethod
    def _files2run(yaml_obj, variables):
        tmpl = yaml_obj.get("files", None)
        if tmpl is None:
            return []
        pattern = jinja2.Template(tmpl).render(variables)
        return sorted(glob.glob(pattern))


class PlayYaml(object):
    def __init__(self, index, files, samples, cmd):
        self._index = index
        self._cmd = cmd
        self._files = files
        self._samples = samples

    def commands(self):
        if len(self._samples) > 0 or len(self._files) > 0:
            for i in range(max(len(self._samples), len(self._files))):
                yield self._cmd.compose(self._index, i)
        else:
            yield self._cmd.compose(self._index, -1)

    def play_run(self, run_index):
        file = self._files[run_index] if len(self._files) > 0  else None
        sample = self._samples[run_index] if len(self._samples) > 0 else None
        return PlayRun(file=file, sample=sample, index=run_index)


class PlayRun(object):
    def __init__(self, file, sample, index):
        self._file = file
        self._sample = sample
        self._index = index

    def vars(self, variables):
        v = copy.deepcopy(variables)
        if self._sample:
            v["curr_sample"] = self._sample
        if self._file:
            v["curr_file"] = self._file
        return v

    def yaml(self, yaml_obj):
        y = copy.deepcopy(yaml_obj)
        y.pop("samples")
        y.pop("files")

    def name(self):
        if self._sample:
            return "sample_" + self._index
        if self._file:
            return os.path.basename(self._file)
        return "None"


class JobCommand(object):
    def __init__(self, playbook_path, tsv_path, var_files, log_dir):
        self._playbook_path = playbook_path
        self._tsv_path = tsv_path
        self._var_files = var_files
        self._log_dir = log_dir

    @staticmethod
    def _next_id(play_index, run_index):
        return job_id_generator.get_next(["play_" + str(play_index) + (str(run_index) if run_index >= 0 else "")])

    def _options(self):
        options = []
        if self._log_dir is not None:
            options.append("--log_dir %s" % self._log_dir)
        if self._tsv_path is not None:
            options.append("--samples %s" % self._tsv_path)
        for var_file in self._var_files:
            options.append("--vars %s" % var_file)

    def compose(self, play_index, run_index):
        executable = "ngspyeasy_play_run"
        cmd = [executable,
               self._playbook_path,
               "--play_index", str(play_index)]
        if run_index >= 0:
            cmd += ["--run_index", str(run_index)]
        cmd += self._options()
        return self._next_id(play_index, run_index), cmd


def parse(playbook_path, tsv_path, var_files, log_dir):
    plays = _read_plays(playbook_path)
    logger().info("Number of plays: %s" % len(plays))

    samples = _read_samples(tsv_path)
    logger().info("Number of samples: %s" % len(samples))

    vars = _read_variables(var_files)
    vars["all_samples"] = samples

    return PlayBookYaml(JobCommand(playbook_path, tsv_path, var_files, log_dir), plays, samples, vars)


def _read_plays(playbook_path):
    logger().info("Reading playbook yaml...")
    with open(playbook_path, 'r') as stream:
        return yaml.load(stream)


def _read_samples(tsv_path):
    if tsv_path is None:
        return []

    logger().debug("TSV config path: %s" % tsv_path)
    tsv_conf = tsv_config.parse(tsv_path)

    logger().info("TSV config first line: %s" % str(tsv_conf.row_at(0)))
    return list(tsv_conf.all_rows())


def _read_variables(var_files):
    d = dict()
    for var_file in var_files:
        with open(var_file, 'r') as stream:
            more_vars = yaml.load(stream)
            d.update(more_vars)
    return d
