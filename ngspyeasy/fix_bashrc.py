#!/usr/bin/env python

import os
import docker


def path_of_fixed_bashrc():
    bashrc = os.path.join(docker.HOME, ".bashrc")
    bashrc_fixed = os.path.join(docker.HOME, ".bashrc_fixed")
    if not os.path.isfile(bashrc):
        return bashrc_fixed

    with open(bashrc, 'r') as input:
        with open(bashrc_fixed, 'a') as output:
            lines = input.readlines()
            ignore = False
            for line in lines:
                output.write('#' + line if ignore else line)
                if line.startswith('# If not running interactively, don\'t do anything'):
                    ignore = True

                if line.startswith('esac') and ignore:
                    ignore = False
    return bashrc_fixed


if __name__ == '__main__':
    os.environ["BASHRC_FIXED"] = path_of_fixed_bashrc()
