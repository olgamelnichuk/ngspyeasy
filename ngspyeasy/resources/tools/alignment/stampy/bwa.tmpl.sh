#!/usr/bin/env bash

time /usr/local/bin/bwa mem \
-M \
-t ${NCPU} \
-R '@RG\tID:${BAM_PREFIX}\tSM:${BAM_PREFIX}\tPU:${PLATFORM_UNIT}\tPL:${NGS_PLATFORM}\tLB:${DNA_PREP_LIBRARY_ID}\tDT:${RUNDATE}' \
${REFFASTA} \
${FQ1} \
${FQ2} | \
samblaster --addMateTags --excludeDups | \
sambamba view -t ${NCPU} -S -f bam /dev/stdin | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${TMP_BAM} /dev/stdin && \
sambamba index ${TMP_BAM}