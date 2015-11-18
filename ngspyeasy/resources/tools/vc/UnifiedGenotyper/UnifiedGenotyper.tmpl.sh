#!/usr/bin/env bash

time java -Xmx6g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R ${REFFASTA} \
-nct ${NCPU} \
-I ${VC_FILTERED_BAM} \
-o ${RAW_VCF} \
-stand_call_conf 30 \
-stand_emit_conf 10 \
--output_mode EMIT_VARIANTS_ONLY \
--dbsnp ${KNOWN_SNPS_B138} \
-dcov 250 \
--unsafe ALL \
--genotype_likelihoods_model BOTH \
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