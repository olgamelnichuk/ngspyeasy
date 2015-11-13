#!/usr/bin/env bash

time /usr/local/bin/bowtie2 ${FAST} --threads ${NCPU} -x ${GENOMEINDEX} \
-I 50 -X 10000 \
--rg-id ${BAM_PREFIX} \
--rg SM:${BAM_PREFIX} \
--rg PU:${NGS_PLATFORM} \
--rg LB:${DNA_PREP_LIBRARY_ID} \
--rg DT:${RUNDATE} \
-1 ${FQ1} -2 ${FQ2} | \
samblaster --addMateTags --excludeDups \
--discordantFile ${DISCORDANT_SAM} \
--splitterFile ${SPLITREAD_SAM} \
--unmappedFile ${UNMAPPED_FASTQ} | \
sambamba view -t ${NCPU} -S -f bam /dev/stdin | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${DUPEMARK_BAM} /dev/stdin && \
sambamba index ${DUPEMARK_BAM} && \
sambamba flagstat -t ${NCPU} ${DUPEMARK_BAM} > ${DUPEMARK_FLAGSTAT_REPORT} && \
bedtools bamtobed -i ${DUPEMARK_BAM} | bedtools merge > ${DUPEMARK_BED_REPORT} && \
sambamba view -t ${NCPU} -S -f bam ${DISCONDANT_SAM} | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${DISCORDANT_BAM} /dev/stdin && \
sambamba index ${DISCORDANT_BAM} && \
sambamba view -t ${NCPU} -S -f bam ${SPLITREAD_SAM} | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${SPLITREAD_BAM} /dev/stdin && \
sambamba index ${SPLITREAD_BAM} && \
rm -v ${DISCORDANT_SAM} && \
rm -v ${SPLITREAD_SAM} && \
rm -rf ${TMP_DIR}/*