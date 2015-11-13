#!/usr/bin/env bash

time python /usr/local/pipeline/Platypus/bin/Platypus.py callVariants \
--nCPU ${NCPU} \
--bamFiles=${FILTERED_BAM} \
--refFile=${REFFASTA} \
--output=${RAW_VCF} \
--logFileName=${RAW_VCF}.log \
--filterDuplicates=${FILTER_DUPLICATES} && \
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