#!/usr/bin/env bash

function transpose() {
    awk '
  {
      for (i=1; i<=NF; i++)  {
          a[NR,i] = $$i
      }
  }
  NF>p { p = NF }
  END {
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' $${1}
}

CHROMS=$$(cut -f 1 ${DUPEMARK_BED} | uniq | transpose)

sambamba view -t ${NCPU} -f bam -F "not (unmapped or mate_is_unmapped) and mapping_quality >=0" ${DUPEMARK_BAM} $${CHROMS} | \
glia -Rru -w 1000 -S 100 -Q 100 -G 4 -f ${REFFASTA} -v ${KNOWN_INDELS} | \
sambamba sort -t ${NCPU} -m 2GB --tmpdir=${TMP_DIR} -o ${DUPEMARK_REALN_BAM} /dev/stdin && \
sambamba index ${DUPEMARK_REALN_BAM} && \
sambamba flagstat -t ${NCPU} ${DUPEMARK_REALN_BAM} > ${DUPEMARK_REALN_FLAGSTAT} && \
bedtools bamtobed -i ${DUPEMARK_REALN_BAM} | bedtools merge > ${DUPEMARK_REALN_BED} && \
rm -rf ${TMP_DIR}/*