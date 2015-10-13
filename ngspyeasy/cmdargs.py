#!/usr/bin/env python

import argparse

import os.path


def existed_directory_path(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError('%s is not an existed directory path' % path)
    return path


def path_basename(path):
    return os.path.basename(path)
