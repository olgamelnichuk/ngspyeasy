#!/usr/bin/env bash

time java -Xmx6g -Djava.io.tmpdir=${TMP_DIR} -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${REFFASTA} \
-nct ${NCPU} \
-I ${VC_FILTERED_BAM} \
-o ${HAPLOTYPE_CALLER_RAW_VCF} \
-stand_call_conf 30 \
-stand_emit_conf 10 \
--output_mode EMIT_VARIANTS_ONLY \
--dbsnp ${KNOWN_SNPS_B138} \
-dcov 250 \
-minPruning 10 \
--unsafe ALL \
-pairHMM VECTOR_LOGLESS_CACHING \
--genotyping_mode DISCOVERY \
--output_mode ${GTMODEGATK} \
--annotation AlleleBalance \
--annotation BaseCounts \
--annotation BaseQualityRankSumTest \
--annotation ChromosomeCounts \
--annotation ClippingRankSumTest \
--annotation Coverage \
--annotation FisherStrand \
--annotation GCContent \
--annotation HaplotypeScore \
--annotation HomopolymerRun \
--annotation InbreedingCoeff \
--annotation LikelihoodRankSumTest \
--annotation LowMQ \
--annotation MappingQualityRankSumTest \
--annotation MappingQualityZero \
--annotation QualByDepth \
--annotation RMSMappingQuality \
--annotation ReadPosRankSumTest \
--annotation SpanningDeletions \
--annotation TandemRepeatAnnotator \
--annotation VariantType && \
time cat ${HAPLOTYPE_CALLER_RAW_VCF} | \
vcffilter -f 'QUAL > 5' -s | \
fix_ambiguous | \
vcfallelicprimitives --keep-geno | \
vcffixup - | \
vcfstreamsort | \
vt normalize -r ${REFFASTA} -q - 2> /dev/null | \
vcfuniqalleles | \
bgzip -c > ${HAPLOTYPE_CALLER_VCF_GZ} && \
tabix ${HAPLOTYPE_CALLER_VCF_GZ} && \
bgzip ${HAPLOTYPE_CALLER_RAW_VCF} && \
tabix ${HAPLOTYPE_CALLER_RAW_VCF_GZ}
