#!/usr/bin/env bash

time /usr/local/pipeline/freebayes/scripts/freebayes-parallel <(/usr/local/pipeline/freebayes/scripts/fasta_generate_regions.py ${REFFASTA} 100000) ${NCPU} \
-f ${REFFASTA} \
-b ${FILTERED_BAM} \
--min-repeat-entropy 1 \
--genotype-qualities > ${RAW_VCF} && \
time cat ${SOUTDocker}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf | \
vcffilter -f 'QUAL > 5' -s | \
fix_ambiguous | \
vcfallelicprimitives --keep-geno | \
vcffixup - | \
vcfstreamsort | \
vt normalize -r ${REFFASTA} -q - 2> /dev/null | \
vcfuniqalleles | \
bgzip -c > ${VCF_GZ} && \
tabix ${VCF_GZ} && \
bgzip  ${RAW_VCF} &&\
tabix ${RAW_VCF_GZ}