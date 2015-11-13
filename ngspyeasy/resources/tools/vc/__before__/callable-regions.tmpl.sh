#!/usr/bin/env bash

time sambamba view \
-t ${NCPU} -F "mapping_quality>=10 and (paired or proper_pair) and (not unmapped)" \
-f bam ${BAM_FILE} | \
bedtools bamtobed | \
bedtools merge > ${MAPPED_READS_BED} && \
time sambamba view \
-t ${NCPU} -F "mapping_quality>1 and (paired or proper_pair) and (not unmapped)" \
-f bam ${BAM_FILE} | \
bedtools genomecov -split -ibam stdin -bga -max 500 | awk '$$4>1' > ${GENOMECOV_BED}