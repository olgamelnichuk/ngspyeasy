# NGSpyeasy
Yet another one NGS pipeline runner based on Ansible automation software.

The basic usage is:

1. Write your pipeline as Ansible playbook.

2. Run it with NGSpyeasy runner

```
$ ngspyeasy --samples /path/to/samples.config.tsv --vars pipeline_vars.yml pipeline_playbook.yml --log_dir /path/to/log_dir
```

where `samples.config.tsv` contains a list of sample configurations your pipeline scripts use:
Here is an example:

```
sample_id   fastq1  fatsq2

sample1 sample1_1.fastq.gz  sample1_2.fastq.gz
sample2 sample2_1.fastq.gz  sample2_2.fastq.gz
sample3 sample3_1.fastq.gz  sample3_2.fastq.gz
...
```
