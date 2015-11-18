#!/usr/bin/env bash

time java -Xmx12g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-nct ${NCPU} \
-R ${REFFASTA} \
-l INFO \
--unsafe ALL \
--validation_strictness SILENT \
-knownSites ${KNOWN_INDELS} \
-knownSites ${KNOWN_SNPS_b138} \
-knownSites ${KNOWN_SNPS_OMNI} \
-knownSites ${KNOWN_SNPS_1000G} \
-I ${BSQR_BAM_IN} \
-o ${RECAL_DATA_TABLE} && \
time java -Xmx12g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T PrintReads \
-nct ${NCPU} \
-R ${REFFASTA} \
-l INFO \
--unsafe ALL \
--validation_strictness SILENT \
--baq RECALCULATE \
--baqGapOpenPenalty 40 \
--BQSR ${RECAL_DATA_TABLE} \
-I ${BSQR_BAM_IN} \
-o ${BSQR_BAM_OUT} && \
sambamba index ${BSQR_BAM_OUT} \
rm -fr ${TMP_DIR}
