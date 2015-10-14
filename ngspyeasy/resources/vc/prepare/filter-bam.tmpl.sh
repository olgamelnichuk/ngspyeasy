#!/usr/bin/env bash

sambamba view \
-t ${NCPU} \
-F "mapping_quality>=20 and proper_pair" \
-f bam -o ${SOUTDocker}/alignments/${FilteredBAM} \
${SOUTDocker}/alignments/${BAMFILE} && \
sambamba index ${SOUTDocker}/alignments/${FilteredBAM} && \
sambamba flagstat -t ${NCPU} ${SOUTDocker}/alignments/${FilteredBAM} > ${SOUTDocker}/alignments/${FilteredBAM}.flagstat