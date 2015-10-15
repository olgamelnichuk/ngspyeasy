#!/usr/bin/env bash

time python /usr/local/pipeline/Platypus/bin/Platypus.py callVariants \
--nCPU ${NCPU} \
--bamFiles=${FILTERED_BAM} \
--refFile=${REFFASTA} \
--output=${RAW_VCF} \
--logFileName=${RAW_VCF}.log \
--filterDuplicates=1 \
--assemble=1 \
--assembleAll=1 \
--assemblyRegionSize=1500 \
--minReads=2  \
--maxGOF=30 \
--bufferSize=50000 \
--maxReads=10000000 \
--minPosterior=5 \
--minMapQual=10 \
--minBaseQual=10 \
--maxSize 10000 \
--maxVariants=8 \
--mergeClusteredVariants=1 \
--maxHaplotypes=50 \
--filterReadsWithDistantMates=1 \
--hapScoreThreshold 10 \
--scThreshold 0.99 \
--filteredReadsFrac 0.7 \
--rmsmqThreshold 20 \
--qdThreshold 0 \
--abThreshold 0.0001 \
--minVarFreq 0.0 && \
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
