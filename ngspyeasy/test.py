#!/usr/bin/env python

import subprocess
import sys
from cmdline_options import run_command
import os
from string import Template

def log(msg):
    print msg


def main():

    print dict(os.environ)

    cmd = "source ~/.profile; echo -n $CLASSPATH\n"

    proc = subprocess.Popen(["bash", "-i", "-c", cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE)
                            #shell=True)
                            #env=dict(os.environ),
                            #executable="/bin/bash")
    (out, err) = proc.communicate()
    output = out
    retcode = 0
    if err:
        log("Command [[\n%s\n]] failed. See logs for details" % cmd)
        output = err
        retcode = proc.returncode

    log(output)
    t  = Template("${A} and ${B} and $C").safe_substitute(A="a", B="b")
    print t


    sys.exit(retcode)




if __name__ == '__main__':
    main()
