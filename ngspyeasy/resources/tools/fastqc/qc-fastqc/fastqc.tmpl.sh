#!/usr/bin/env bash

/usr/local/pipeline/FastQC/fastqc --threads 2 --extract --dir ${TMP_DIR} --outdir ${FASTQ_DIR} ${FQ1} ${FQ2}