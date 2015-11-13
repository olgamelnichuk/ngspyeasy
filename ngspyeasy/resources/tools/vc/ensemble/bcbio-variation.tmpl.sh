#!/usr/bin/env bash

bcbio-variation-recall ensemble \
-c ${NCPU} \
-n 2 \
--nofiltered \
${VCF_GZ} \
${REFFASTA} \
${FREEBAYERS_VCF} \
${HAPLOTYPE_CALLER_VCF}