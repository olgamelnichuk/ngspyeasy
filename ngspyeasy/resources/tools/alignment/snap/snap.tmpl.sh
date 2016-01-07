#!/usr/bin/env bash

time /usr/local/pipeline/snap/snap-aligner paired \
${REFDIR} \
${TRIMMED_FQ1} ${TRIMMED_FQ2} \
-t ${NCPU} \
-b \
-M \
-s 50 1000 \
-H 300000 \
-h 300 \
-d 15 \
-mcp 6000000 \
-map \
-pre \
-I \
-R '@RG\tID:${BAM_PREFIX}\tSM:${BAM_PREFIX}\tPU:${PLATFORM_UNIT}\tPL:${NGS_PLATFORM}\tLB:${DNA_PREP_LIBRARY_ID}\tDT:${RUNDATE}' \
-o -sam -  | \
samblaster --addMateTags --excludeDups \
--discordantFile ${DISCORDANT_SAM} \
--splitterFile ${SPLITREAD_SAM} \
--unmappedFile ${UNMAPPED_FASTQ} | \
sambamba view -t ${NCPU} -S -f bam /dev/stdin | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${DUPEMARK_BAM} /dev/stdin && \
sambamba index ${DUPEMARK_BAM} && \
sambamba flagstat -t ${NCPU} ${DUPEMARK_BAM} > ${DUPEMARK_FLAGSTAT} && \
bedtools bamtobed -i ${DUPEMARK_BAM}| bedtools merge > ${DUPEMARK_BED} && \
sambamba view -t ${NCPU} -S -f bam ${DISCORDANT_SAM} | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${DISCORDANT_BAM} /dev/stdin && \
sambamba index ${DISCORDANT_BAM} && \
sambamba view -t ${NCPU} -S -f bam ${SPLITREAD_SAM} | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${SPLITREAD_BAM} /dev/stdin && \
sambamba index ${SPLITREAD_BAM} && \
rm -v ${DISCORDANT_SAM} && \
rm -v ${SPLITREAD_SAM} && \
rm -rf ${TMP_DIR}/*