#!/usr/bin/env bash

/usr/local/pipeline/FastQC/fastqc --threads 4 \
--extract --quiet --dir ${TMP_DIR} --outdir ${FASTQ_DIR} ${PAIRED_FQ1} ${PAIRED_FQ2} ${UNPAIRED_FQ1} ${UNPAIRED_FQ2}