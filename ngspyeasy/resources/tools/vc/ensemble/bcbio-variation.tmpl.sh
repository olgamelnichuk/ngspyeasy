#!/usr/bin/env bash

bcbio-variation-recall ensemble \
-c ${NCPU} \
-n 2 \
--nofiltered \
${VCF_GZ} \
${REFFASTA} \
${FREEBAYERS_VCF_GZ} \
${HAPLOTYPE_CALLER_VCF_GZ}