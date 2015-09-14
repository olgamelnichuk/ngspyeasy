#!/usr/bin/env python

import os.path

def get_resource_path(resource_name):
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, "resources", resource_name)

