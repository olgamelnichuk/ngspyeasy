#!/usr/bin/env bash

java -XX:ParallelGCThreads=1 -jar /usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
-threads ${NCPU} \
${FQ1} ${FQ2} \
${PAIRED_FQ1} ${UNPAIRED_FQ1} \
${PAIRED_FQ2} ${UNPAIRED_FQ2} \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
AVGQUAL:2 \
MINLEN:75