#!/usr/bin/env bash

time java -Xmx12g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/pipeline/picardtools/${PICARD_VERSION}/picard.jar CleanSam \
TMP_DIR=${TMP_DIR} \
CREATE_INDEX=FALSE \
VALIDATION_STRINGENCY=SILENT \
INPUT=${STAMPY_DUPEMARK_TMP_BAM} \
OUTPUT=${STAMPY_DUPEMARK_CLEANSAM_TMP_BAM} && \
sambamba index ${STAMPY_DUPEMARK_CLEANSAM_TMP_BAM}