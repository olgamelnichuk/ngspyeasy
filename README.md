# NGSPyEasy
Python version of [NGSeasy](https://github.com/KHP-Informatics/ngseasy)

The basic usage is: 
```
$ python /<path-to-ngspyeasy>/ngspyeasy.py -c ngseasy_test.config.tsv -d /home/bob/ngs_projects
```
where config file is compartible with the one you use in NGSeasy, but NGSPyEasy ignores `PROJECT_DIR` column and uses `ngs_projects` directory instead (just for simplicity). 

`ngspyeasy` doesn't require a specific directory to change to before the start; config file option is just a name of the config file in the `projects_home/config_files` directory, but it accepts the full path as well.

The long command line option names are also supported:
```
$ python /opt/ngspyeasy/bin/ngspyeasy.py --config ngseasy_test.config.tsv \
    --projects-dir /home/bob/ngs_projects \
    --resources-dir /nfs/ngs_resources --verbose
```

Running just one step is also possible. Here is an example of running `trimmomatic` step:
```
$ python /opt/ngspyeasy/bin/ngspyeasy.py trimmomatic --config ngseasy_test.config.tsv \
    --projects-dir /home/bob/ngs_projects \
    --resources-dir /nfs/ngs_resources --verbose
```

