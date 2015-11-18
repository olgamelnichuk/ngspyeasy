#!/usr/bin/env bash

time java -Xmx6g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${REFFASTA} \
-nct ${NCPU} \
-I ${VC_FILTERED_BAM} \
-o ${RAW_VCF} \
--dbsnp ${KNOWN_SNPS_B138} \
--unsafe ALL \
-stand_call_conf 30 \
-stand_emit_conf 10 \
--output_mode EMIT_VARIANTS_ONLY \
-dcov 250 \
-minPruning 10 \
--unsafe ALL \
-pairHMM VECTOR_LOGLESS_CACHING \
--genotyping_mode DISCOVERY && \
time cat ${RAW_VCF} | \
vcffilter -f 'QUAL > 5' -s | \
fix_ambiguous | \
vcfallelicprimitives --keep-geno | \
vcffixup - | \
vcfstreamsort | \
vt normalize -r ${REFFASTA} -q - 2> /dev/null | \
vcfuniqalleles | \
bgzip -c > ${VCF_GZ} && \
tabix ${VCF_GZ} && \
bgzip ${RAW_VCF} && \
tabix ${RAW_VCF_GZ}
