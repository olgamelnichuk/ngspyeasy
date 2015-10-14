#!/usr/bin/env bash

time sambamba view \
-t ${NCPU} -F "mapping_quality>=10 and (paired or proper_pair) and (not unmapped)" \
-f bam ${SOUTDocker}/alignments/${BAMFILE} | \
bedtools bamtobed | \
bedtools merge > ${SOUTDocker}/reports/${BAMFILE}.mapped.reads.bed && \
time sambamba view \
-t ${NCPU} -F "mapping_quality>1 and (paired or proper_pair) and (not unmapped)" \
-f bam ${SOUTDocker}/alignments/${BAMFILE} | \
bedtools genomecov -split -ibam stdin -bga -max 500 | awk '\$4>1' > ${SOUTDocker}/reports/${BAMFILE}.genomecov.bed

# quick fix to get mapped regrions file for freebayes-parallel
# START Make regions files for var-callers
awk '{print $1":"$2".."$3}' ${SOUT}/reports/${BAMFILE}.mapped.reads.bed > ${SOUT}/reports/${BAMFILE}.mapped.regions
awk '{print $1":"$2"-"$3}'  ${SOUT}/reports/${BAMFILE}.mapped.reads.bed > ${SOUT}/reports/${BAMFILE}.platypus.intervals