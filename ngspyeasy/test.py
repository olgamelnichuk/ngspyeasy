import subprocess
import sys
import os


def log(msg):
    print msg


def main():
    cmd = "source ~/.bashrc && env"

    proc = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True,
                            env=os.environ.copy(),
                            executable="/bin/sh")
    (out, err) = proc.communicate()
    output = out
    retcode = 0
    if err:
        log("Command [[\n%s\n]] failed. See logs for details" % cmd)
        output = err
        retcode = proc.returncode

    log("cmd: \n" + output)
    sys.exit(retcode)


if __name__ == '__main__':
    main()
