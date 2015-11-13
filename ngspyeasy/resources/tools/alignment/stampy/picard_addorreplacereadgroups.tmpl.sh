#!/usr/bin/env bash

time java -Xmx12g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/pipeline/picardtools/${PICARD_VERSION}/picard.jar AddOrReplaceReadGroups \
TMP_DIR=${TMP_DIR} \
CREATE_INDEX=FALSE \
VALIDATION_STRINGENCY=SILENT \
INPUT=${DUPEMARK_TMPCLEANSAM_BAM} \
OUTPUT=${DUPEMARK_BAM}  \
RGID=${BAM_PREFIX}\tSM:${BAM_PREFIX} \
RGLB=${DNA_PREP_LIBRARY_ID} \
RGPL=${NGS_PLATFORM} \
RGPU=${PLATFORM_UNIT} \
RGSM=${BAM_PREFIX} \
RGDT=${RUNDATE} && \
sambamba index ${DUPEMARK_BAM} && \
sambamba flagstat -t ${NCPU} ${DUPEMARK_BAM} > ${DUPEMARK_FLAGSTAT_REPORT} && \
bedtools bamtobed -i ${DUPEMARK_BAM} | bedtools merge > ${DUPEMARK_BED} && \
rm -v ${DUPEMARK_TMPCLEANSAM_BAM} && \
rm -v ${DUPEMARK_TMP_BAM} && \
rm -v ${DUPEMARK_TMP_BAM}.bai && \
rm -v ${TMP_BAM}.bai && \
rm -v ${TMP_BAM} && \
rm -v ${DUPEMARK_TMPCLEANSAM_BAM}.bai