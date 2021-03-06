# annotate segments / peaks

set -e
source `dirname $0`/../config.sh
verbose=0
# input segmentation (both strands with segment IDs)
INSEGM=$1 
nsegm=$(wc -l ${INSEGM} | awk '{print $1}')
if [ "${nsegm}" -eq 0 ]; then
   INSEGMBASE=${INSEGM}
   >${INSEGMBASE}.annot.final
   >${INSEGMBASE}.unannotated.bed
   echo "SKIPPING: ${INSEGM} is empty"
   exit
fi
# col.5 = segment ID
# col.4 = tissue name
#ANNOTDIR=${SPARPATH}/annot
ANNOTDIR=`dirname $0`/../annot
ANNOTS1=${ANNOTDIR}/hsa19.fulltable.no_mRNA_no_lncRNA_no_piRNA_nomiRP.unique_LOC.bed
ANNOTS2=${ANNOTDIR}/hsa19.fulltable.miRPonly.unique_LOC.bed
ANNOTS3=${ANNOTDIR}/hsa19.fulltable.no_mRNA_no_lncRNA_piRNAonly.unique_LOC.bed
ANNOTS4=${ANNOTDIR}/hsa19.fulltable.no_mRNA_no_lncRNA_rRNAonly.unique_LOC.bed
ANNOTS5=${ANNOTDIR}/hsa19.fulltable.no_mRNA_no_lncRNA.unique_LOC.bed
#INSEGMBASE=`basename ${INSEGM}`
INSEGMBASE=${INSEGM}
if [ ${verbose} -eq 1 ]; then echo "${INSEGMBASE}"; fi
# step1.
# 1.1 overlap with miRNA-5p/3p/5p3pno, snoRNA, snRNA, scRNA, tRF3, tRF5, tRNA

${BEDTOOLS} intersect -a ${INSEGM}  -b ${ANNOTS1} -wao -s -f 0.8 > ${INSEGMBASE}.step1.bed

# 1.2 col. 8 = -1 put for the next step
numFields=$(awk '{print NF; exit}' ${INSEGM})

if [ ${verbose} -eq 1 ]; then echo "Number of fields: ${numFields}"; fi
overlapField=$(awk 'BEGIN{print '${numFields}'+7;}')
if [ ${verbose} -eq 1 ]; then echo "Overlap field: ${overlapField}"; fi

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep2.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          overlapField=numInputFields+7;
          annotChrStart=numInputFields+2;
          annotChrEnd=numInputFields+3; 
          $overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step1.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
           awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step1.annot.out

# step 2.
# overlap with miRNA precursor

if [ ! -f "${INSEGMBASE}.forstep2.bed" ]; then
  #exit, if nothing to do in step 2
  echo "SKIPPING: nothing to do for STEP 2"
  exit
fi

${BEDTOOLS} intersect -a ${INSEGMBASE}.forstep2.bed -b ${ANNOTS2} -wao -s -f 0.8 > ${INSEGMBASE}.step2.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep3.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          overlapField=numInputFields+7;
          annotChrStart=numInputFields+2;
          annotChrEnd=numInputFields+3; 
          $overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step2.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr > ${INSEGMBASE}.step2.annot.out 
 
# step 3.
# overlap with piRNA

if [ ! -f "${INSEGMBASE}.forstep3.bed" ]; then
  #exit, if nothing to do in step 3
  echo "SKIPPING: nothing to do for STEP 3"
  exit
fi



${BEDTOOLS} intersect -a ${INSEGMBASE}.forstep3.bed -b ${ANNOTS3} -wao -s -f 0.8 > ${INSEGMBASE}.step3.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep4.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          #overlapField=numInputFields+7;
          #annotChrStart=numInputFields+2;
          #annotChrEnd=numInputFields+3; 
          #$overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step3.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
     awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step3.annot.out

# step 4.
# annotate rRNA

if [ ! -f "${INSEGMBASE}.forstep4.bed" ]; then
  #exit, if nothing to do in step 4
  echo "SKIPPING: nothing to do for STEP 4"
  exit
fi



${BEDTOOLS} intersect -a ${INSEGMBASE}.forstep4.bed -b ${ANNOTS4} -wao -s -f 0.8 > ${INSEGMBASE}.step4.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.forstep5.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          #overlapField=numInputFields+7;
          #annotChrStart=numInputFields+2;
          #annotChrEnd=numInputFields+3; 
          #$overlapField=$overlapField/($annotChrEnd-$annotChrStart);
          print;
       }
     }' ${INSEGMBASE}.step4.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
           awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step4.annot.out


# step 5. 
# re-annotate remaining segments/peaks without restriction on the overlapping criteria

if [ ! -f "${INSEGMBASE}.forstep5.bed" ]; then
  #exit, if nothing to do in step 5
  echo "SKIPPING: nothing to do for STEP 5"
  exit
fi

${BEDTOOLS} intersect -a ${INSEGMBASE}.forstep5.bed -b ${ANNOTS5} -wao -s > ${INSEGMBASE}.step5.bed 

awk 'BEGIN{OFS="\t"; numInputFields='${numFields}'+0;annotStartIdx=numInputFields+2;nextStepFile="'${INSEGMBASE}'.unannotated.bed"}{
       
       if ($annotStartIdx=="-1") {for (i=1; i<numInputFields;i++) printf "%s\t", $i > nextStepFile; printf "%s\n", $numInputFields > nextStepFile;}
       else
       {
          # overlap percentage
          overlapField=numInputFields+7;
          inputChrStart=2;
          inputChrEnd=3; 
          $overlapField=$overlapField/($inputChrEnd-$inputChrStart);
          print;
       }
     }' ${INSEGMBASE}.step5.bed | sort -k1,1 -k2,2n -k3,3n -k10,10 -k${overlapField},${overlapField}nr | \
           awk 'BEGIN{OFS="\t"}{
                  segmID=$5;
                  if (prev_segmID != segmID)
                  {
                     print;
                  }
                  prev_segmID = $5;
                }' > ${INSEGMBASE}.step5.annot.out

#cat ${INSEGMBASE}.*.annot.out > ${INSEGMBASE}.annot.final
cat ${INSEGMBASE}.*.annot.out | sort -k1,1 -k2,2n -k3,3n -k6,6 > ${INSEGMBASE}.annot.final
#rm ${INSEGMBASE}*step*.*
