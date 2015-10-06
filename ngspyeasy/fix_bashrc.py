#!/usr/bin/env python

def main():
    with open('~/.bashrc', 'r') as input:
        with open('~/.bashrc_fixed', 'a') as output:
            lines = input.readlines()
            ignore = False
            for line in lines:
                output.write('#' + line if ignore else line)
                if line.startswith('# If not running interactively, don\'t do anything'):
                    ignore = True

                if line.startswith('esac') and ignore:
                    ignore = False


if __name__ == '__main__':
    main()
