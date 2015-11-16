#!/usr/bin/env bash

sambamba view \
-t ${NCPU} \
-F "mapping_quality>=20 and proper_pair" \
-f bam -o ${VC_FILTERED_BAM} \
${VC_BAM_IN} && \
sambamba index ${VC_FILTERED_BAM} && \
sambamba flagstat -t ${NCPU} ${VC_FILTERED_BAM} > ${VC_FILTERED_BAM}.flagstat