#!/usr/bin/env bash

wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip

unzip /tmp/Trimmomatic-0.32.zip -d /tmp

java -XX:ParallelGCThreads=1 -jar /tmp/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
-threads ${NCPU} \
${FQ1} ${FQ2} \
${PAIRED_FQ1} ${UNPAIRED_FQ1} \
${PAIRED_FQ2} ${UNPAIRED_FQ2} \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
AVGQUAL:2 \
MINLEN:75 \
ILLUMINACLIP:${ADAPTER_FA}:2:30:10:5:true