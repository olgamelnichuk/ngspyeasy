#!/usr/bin/env bash

sambamba view \
-t ${NCPU} \
-F "mapping_quality>=20 and proper_pair" \
-f bam -o ${FILTERED_BAM} \
${BAM_FILE} && \
sambamba index ${FILTERED_BAM} && \
sambamba flagstat -t ${NCPU} ${FILTERED_BAM} > ${FILTERED_BAM}.flagstat