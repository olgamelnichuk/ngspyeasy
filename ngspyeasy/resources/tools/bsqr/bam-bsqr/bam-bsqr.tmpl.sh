#!/usr/bin/env bash

time bam recab \
--in ${BSQR_BAM_IN} \
--out -.bam \
--refFile ${REFFASTA} \
--dbsnp ${DBSNP_RECAB} \
--storeQualTag OQ \
--maxBaseQual 40 | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${BAM_OUT} /dev/stdin && \
sambamba index ${BSQR_BAM_OUT} && \
rm -rf ${TMP_DIR}
