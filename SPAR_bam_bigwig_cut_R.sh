

source `dirname $0`/config.sh

set -e

INFILE=$1

if [ $# -lt 1 ]
then
  echo "USAGE: `basename $0` <reads.fastq|aln.bam|signal.bigWig> <SPAR_output_directory>"
  exit 1
fi


if [ ! -f "${INFILE}" ]; then
  echo -e "*****ERROR: FASTQ file\n${INFILE}\nnot found!"
  exit 1
fi

# determine file type
# if BAM, skip mapping, proceed
# if FASTQ, map, then proceed from BAM to BEDGRAPH conversion
# if BIGWIG, proceed from segmentation

FILETYPE=${INFILE##*.}
lcFILETYPE=$(echo ${FILETYPE} | awk '{print tolower($0)}') # lower-case 
ucFILETYPE=$(echo ${FILETYPE} | awk '{print toupper($0)}') # upper-case 

isBAM=0
isBIGWIG=0
isFASTQ=0
#if [ "${FILETYPE,,}" = "bam" ]; then
if [ "${lcFILETYPE}" = "bam" ]; then
  isBAM=1
fi

#if [ "${FILETYPE,,}" = "bigwig" ]; then # this is for bash 4+ only
if [ "${lcFILETYPE}" = "bigwig" ]; then
  isBIGWIG=1
fi

if [ "${lcFILETYPE}" = "fastq" ]; then
  isFASTQ=1
fi


EXPNAME=`basename "${INFILE}" ".${FILETYPE}"`
#EXPNAME=${EXPNAME%_*}


# set SPARDIR output directory to specified directory if any, or use
# default directory from config file
 
SPARDIR=${2:-${SPARDIR}} 
SPARPATH=`dirname $0`


maxReadLength=${3:-44}

if [ "${isFASTQ}" = 1 ]; then
  OUTDIR=${SPARDIR}/${EXPNAME}_m${maxMismatchCnt}_map${maxMapCnt}
  OUTBAM=${OUTDIR}/Aligned.out.filtered.hardClipped.sorted.bam
elif [ "${isBAM}" = 1 ]; then 
  # if input is mapped BAM file
  if [ -z "$2" ]; then
    # if no output directory is provided, use BAM file directory
    OUTDIR=`dirname ${INFILE}`  #${SPARDIR}/${EXPNAME}
  else
    # if output directory is specified, use it
    OUTDIR=$2
  fi
  OUTBAM=${INFILE}
elif [ "${isBIGWIG}" = 1 ]; then
  if [ -z "$2" ]; then
    # if no output directory is provided, use bigWig file directory
    OUTDIR=`dirname ${INFILE}`  #${SPARDIR}/${EXPNAME}
  else
    # if output directory is specified, use it
    OUTDIR=$2
  fi
  OUTBAM=${INFILE}
fi

echo "FILETYPE=${FILETYPE}; isBAM=${isBAM}; EXPNAME=${EXPNAME}; OUTDIR=${OUTDIR}; OUTBAM=${OUTBAM}"
mkdir -p ${OUTDIR}

#LOGSPAR=${OUTDIR}/SPAR.log
#LOGSPAR=${OUTDIR}/`basename ${INFILE}`.SPAR.log
LOGSPAR=${OUTDIR}/`basename ${OUTBAM}`.SPAR.log
>${LOGSPAR}

function runScript
{
  bash ${SPARPATH}/scripts/$1
}

function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1" | tee -a ${LOGSPAR}
}

function printL
{
  echo -e "$1" | tee -a ${LOGSPAR}
}

printL "Output directory:\n${OUTDIR}\n"

printT "SPAR ${ucFILETYPE} run started"

if [ "${isFASTQ}" = 1 ]; then
  runScript "run_star_smrna2.sh ${INFILE} ${maxMismatchCnt} ${maxMapCnt} ${OUTDIR}"
fi

if [ "${isBAM}" = 1 ] || [ "${isFASTQ}" = 1 ]; then
  printT "Converting BAM to bedGraph"
  runScript "bam_to_bedgraph2_cut.sh ${OUTBAM} ${OUTDIR} ${maxReadLength}"

  if [ "${isBAM}" = 1 ]; then
    OUTBAM=${OUTDIR}/`basename "${INFILE}"`
  fi
  #echo "OUTBAM=${OUTBAM}"
  printT "Converting bedGraph to bigWig"
  runScript "bedgraph_to_bigwig.sh ${OUTBAM}.pos.bedgraph"
  runScript "bedgraph_to_bigwig.sh ${OUTBAM}.neg.bedgraph"
fi

if [ "${isBIGWIG}" = 1 ]; then
  OUTBAM=${OUTDIR}/`basename "${INFILE}"`
  # convert bigwig to bedGraph
  INBIGWIG=${INFILE}
  # need to SET the strand for BIGWIG!!!!
  strand=${INBIGWIG%.*}
  strand=${strand##*.}
  if [ "${strand}" != "pos" ] && [ "${strand}" != "neg" ]; then
     # set strand="positive" if strand information is not provided
     strand="pos"
  fi
  #echo "strand=${strand}"
  OUTBEDGRAPH=${OUTDIR}/`basename "${INBIGWIG}"`.${strand}.bedgraph
  ${BIGWIGTOBEDGRAPH} ${INBIGWIG} ${OUTBEDGRAPH} 
  printT "Segmenting BIGWIG"
  runScript "segment_bedgraph_entropy_position2.sh ${OUTBEDGRAPH} ${strand}"
  numFields=$(awk '{print NF; exit}' ${OUTBEDGRAPH}.segm)
  printT "Creating called peaks track BIGWIG"
  runScript "create_peak_track.sh ${OUTBEDGRAPH}.segm"
  printT "Annotating peaks BIGWIG"
  runScript "annotate_segm3.sh ${OUTBEDGRAPH}.segm"
else

printT "Segmenting [positive strand]"
#runScript "segment_bedgraph_entropy.sh ${OUTBAM}.pos.bedgraph pos"
runScript "segment_bedgraph_entropy_position2.sh ${OUTBAM}.pos.bedgraph pos"
printT "Segmenting [negative strand]"
#runScript "segment_bedgraph_entropy.sh ${OUTBAM}.neg.bedgraph neg"
runScript "segment_bedgraph_entropy_position2.sh ${OUTBAM}.neg.bedgraph neg"


numFields=$(awk '{print NF; exit}' ${OUTBAM}.pos.bedgraph.segm)

printT "Creating track files (Raw signal)"
runScript "create_bedgraph_track.sh ${OUTBAM}.pos.bedgraph"
runScript "create_bedgraph_track.sh ${OUTBAM}.neg.bedgraph"

printT "Creating track files (Called peaks)"
runScript "create_peak_track.sh ${OUTBAM}.pos.bedgraph.segm"
runScript "create_peak_track.sh ${OUTBAM}.neg.bedgraph.segm"


printT "Annotating [positive strand]"
runScript "annotate_segm3.sh ${OUTBAM}.pos.bedgraph.segm"
printT "Annotating [negative strand]"
runScript "annotate_segm3.sh ${OUTBAM}.neg.bedgraph.segm"


fi


printT "DONE."

#chr,a,b,peakID,exprVal,strand,entropy5p,normEntropy5p,maxEntropy5p, entropy3p, normEntropy3p, maxEntropy3p, maxProb5p, maxProb3p;

# add headers
#awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction", "max5pPosition", "max3pPosition", "annotChr", "annotChrStart", "annotChrEnd", "annotID", "annotRNAclass","annotStrand","overlap"; exit}' |  cat - ${OUTBAM}.*.bedgraph.segm.annot.final > ${OUTBAM}.annot.final


#awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction", "max5pPosition", "max3pPosition", "annotChr", "annotChrStart", "annotChrEnd", "annotID", "annotRNAclass","annotStrand","overlap"; exit}' > ${OUTBAM}.annot.final
#chr,a,b,peakID,exprVal,strand,entropy5p,normEntropy5p,maxEntropy5p, entropy3p, normEntropy3p, maxEntropy3p, maxProb5p, maxProb3p, maxPosition5p, maxPosition3p, beforePeakReadCnt, peakFoldChange;
awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction", "max5pPosition", "max3pPosition", "beforePeakReadCnt", "peakFoldChange", "annotChr", "annotChrStart", "annotChrEnd", "annotID", "annotRNAclass","annotStrand","overlap"; exit}' > ${OUTBAM}.annot.final
positiveannot="${OUTBAM}.pos.bedgraph.segm.annot.final"
if [ -f ${positiveannot} ]; then
  cat ${positiveannot} >> ${OUTBAM}.annot.final
fi
negativeannot="${OUTBAM}.neg.bedgraph.segm.annot.final"
if [ -f "${negativeannot}" ]; then
  cat ${negativeannot} >> ${OUTBAM}.annot.final
fi


#awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction", "max5pPosition", "max3pPosition"; exit}' |  cat - ${OUTBAM}.*.bedgraph.segm.unannotated.bed > ${OUTBAM}.unannot
#awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction", "max5pPosition", "max3pPosition"; exit}' > ${OUTBAM}.unannot
awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction", "max5pPosition", "max3pPosition", "beforePeakReadCnt", "peakFoldChange"; exit}' > ${OUTBAM}.unannot
positiveunannot="${OUTBAM}.pos.bedgraph.segm.unannotated.bed"
if [ -f "${positiveunannot}"  ]; then
  cat ${positiveunannot} >> ${OUTBAM}.unannot
fi

negativeunannot="${OUTBAM}.neg.bedgraph.segm.unannotated.bed"
if [ -f "${negativeunannot}" ]; then
  cat ${negativeunannot} >> ${OUTBAM}.unannot
fi


#cat ${OUTBAM}.*.bedgraph.segm.annot.final > ${OUTBAM}.annot.final
#cat ${OUTBAM}.*.bedgraph.segm.unannotated.bed > ${OUTBAM}.unannot

#echo "numFields=${numFields}"

# annotation summary
annotSummary=${OUTBAM}.mapped_reads_annotation_summary.txt
#awk 'BEGIN{OFS="\t";}{if ($0~/^#/) {next}; if (NR==FNR) {exprVal=$5; rnaClass=$21; exprPerClass[rnaClass]+=exprVal;classCnt[rnaClass]+=1;totalExprAnnot+=exprVal; totalAnnotPeakCnt+=1}else{exprVal=$5; totalExprUnannot+=exprVal;totalUnannotPeakCnt+=1;}}END{totalPeakCnt=totalAnnotPeakCnt+totalUnannotPeakCnt; totalExpr=totalExprAnnot+totalExprUnannot; for (rnaClass in exprPerClass) print rnaClass, classCnt[rnaClass], exprPerClass[rnaClass], exprPerClass[rnaClass]/totalExpr; print "Annotated", totalAnnotPeakCnt, totalExprAnnot, totalExprAnnot/totalExpr; print "Unannotated",totalUnannotPeakCnt,totalExprUnannot,totalExprUnannot/totalExpr}' ${OUTBAM}.annot.final ${OUTBAM}.unannot | sort -k1,1 | awk 'BEGIN{OFS="\t"; print "#RNA class","Peaks","Reads","Fraction of reads"}{print}' > ${annotSummary} 
awk 'BEGIN{OFS="\t";totalAnnotPeakCnt=0; totalUnannotPeakCnt=0;totalExprAnnot=0;totalExprUnannot=0;numFields='${numFields}'+0;}{if ($0~/^#/) {next}; if (NR==FNR) {exprVal=$5; rnaClass=$(numFields+5); exprPerClass[rnaClass]+=exprVal;classCnt[rnaClass]+=1;totalExprAnnot+=exprVal; totalAnnotPeakCnt+=1}else{exprVal=$5; totalExprUnannot+=exprVal;totalUnannotPeakCnt+=1;}}END{totalPeakCnt=totalAnnotPeakCnt+totalUnannotPeakCnt; totalExpr=totalExprAnnot+totalExprUnannot; for (rnaClass in exprPerClass) print rnaClass, classCnt[rnaClass], exprPerClass[rnaClass], exprPerClass[rnaClass]/totalExpr; propAnnot=0; if (totalExpr>0) propAnnot=totalExprAnnot/totalExpr; print "Annotated", totalAnnotPeakCnt, totalExprAnnot, propAnnot; propUnannot=0; if (totalExpr>0) propUnannot=totalExprUnannot/totalExpr; print "Unannotated",totalUnannotPeakCnt,totalExprUnannot,propUnannot}' ${OUTBAM}.annot.final ${OUTBAM}.unannot | sort -k1,1 | awk 'BEGIN{OFS="\t"; print "#RNA class","Peaks","Reads","Fraction of reads"}{print}' > ${annotSummary} 
finalAnnot="${OUTBAM}.annot.final"
finalUnannot="${OUTBAM}.unannot"

awk 'BEGIN{FS="\t"; OFS="\t";numFields='${numFields}'+0;}
     {
       if (NR==1) next; # skip header
       annotID=$(numFields+4);
       #print annotID;
       n=split(annotID,a,":");
       geneID=a[n];
       annotClass=$(numFields+5);
       if (ids[geneID]!=1)
       {
          geneCnt[annotClass]+=1;
          ids[geneID]=1;
       }
     }
     END{ for (c in geneCnt)
           print c, geneCnt[c]
     }' ${finalAnnot} | sort -k1,1 | awk 'BEGIN{OFS="\t"; FS="\t"; print "#RNA class", "Genes"}{print}' >> ${annotSummary}

annotLengthSummary=${OUTBAM}.annot.peak.length.stats
>${annotLengthSummary}
#echo -e "\nLength of annotated peaks:\n" >> ${annotLengthSummary}
awk 'BEGIN{FS="\t"; OFS="\t"}
     {
       if (NR==1) next;
       peakStart = $2; # 0-based
       peakEnd = $3; # 1-based according to 0-based, half-open UCSC notation
       peakLength = peakEnd-peakStart;
       lengthCnt[peakLength]+=1;  
       totalPeakCnt+=1;
     }
     END{ for (l in lengthCnt)
          {
            lengthProp = lengthCnt[l]/totalPeakCnt
            print l, lengthCnt[l], lengthProp;
          }
     }' ${finalAnnot} | sort -k1,1n | awk 'BEGIN{FS="\t"; OFS="\t"; print "#PeakLength", "Count", "Fraction"}{print}' > ${annotLengthSummary}

unannotLengthSummary=${OUTBAM}.unannot.peak.length.stats
>${unannotLengthSummary}
#echo -e "\nLength of annotated peaks:\n" >> ${annotLengthSummary}
awk 'BEGIN{FS="\t"; OFS="\t"}
     {
       if (NR==1) next;
       peakStart = $2; # 0-based
       peakEnd = $3; # 1-based according to 0-based, half-open UCSC notation
       peakLength = peakEnd-peakStart;
       lengthCnt[peakLength]+=1;  
       totalPeakCnt+=1;
     }
     END{ for (l in lengthCnt)
          {
            lengthProp = lengthCnt[l]/totalPeakCnt
            print l, lengthCnt[l], lengthProp;
          }
     }' ${finalUnannot} | sort -k1,1n | awk 'BEGIN{FS="\t"; OFS="\t"; print "#PeakLength", "Count", "Fraction"}{print}' > ${unannotLengthSummary}





printL "\n===Output==="
printL "Output directory: ${OUTDIR}"
printL "LOG file: ${LOGSPAR}"
if [ "${isFASTQ}" = 1 ]; then
  printL "\nMapping output:"
  printL "${OUTBAM}"
fi

if [ "${isBAM}" = 1 ] || [ "${isFASTQ}" = 1 ]; then
  printL "\nTrack files (Raw signal):"
  ls ${OUTBAM}*.bigWig | tee -a ${LOGSPAR}
fi
printL "\nTrack files (Called peaks):"
ls ${OUTBAM}*.bigBed | tee -a ${LOGSPAR}

printL "\nAnnotation output:"
#ls ${OUTBAM}.*.bedgraph.segm.annot.final | tee -a ${LOGSPAR}
ls ${OUTBAM}.annot.final | tee -a ${LOGSPAR}

printL "\nUn-annotated output:" 
#ls ${OUTBAM}.*.bedgraph.segm.unannotated.bed | tee -a ${LOGSPAR}
ls ${OUTBAM}.unannot | tee -a ${LOGSPAR}


printL "\nAnnotation summary:"
#ls ${OUTDIR}/mapped_reads_annotation_summary.txt | tee -a ${LOGSPAR}
ls ${annotSummary}  | tee -a ${LOGSPAR}

if [ "${isFASTQ}" = 1 ]; then 
  printL "\nMapping stats:"
  ls ${OUTDIR}/MAPSTAT.txt | tee -a ${LOGSPAR}
  printL "\nMapped reads alignment summary:"
  ls ${OUTDIR}/cigar.stat | tee -a ${LOGSPAR}
fi



printL "\n\n===Run summary==="

if [ "${isFASTQ}" = 1 ]; then 
  printL "FASTQ: ${INFILE}"
  grep -e "Reads \[all\]" ${OUTDIR}/MAPSTAT.txt | awk 'BEGIN{FS="\t"}{printf "Mapped reads: %d [%.4f%%]\n", $2, $3}' | tee -a ${LOGSPAR}
elif [ "${isBAM}" = 1 ]; then
  printL "BAM: ${INFILE}"
elif [ "${isBIGWIG}" = 1 ]; then
  printL "BIGWIG: ${INFILE}"
fi

numAnnot=$(cat ${OUTBAM}.*.bedgraph.segm.annot.final | wc -l)
printL "Annotated loci count: ${numAnnot}"
#printL "Annotated loci by RNA class:"
#cat ${OUTBAM}.*.bedgraph.segm.annot.final | cut -f 19 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$1}' | sort -k2,2nr | tee -a ${LOGSPAR}

numUnannot=$(cat ${OUTBAM}.*.bedgraph.segm.unannotated.bed | wc -l)
printL "Un-annotated loci count: ${numUnannot}"


printL "\nAnnotation summary:"
#ls ${OUTDIR}/mapped_reads_annotation_summary.txt | tee -a ${LOGSPAR}
ls ${annotSummary} | tee -a ${LOGSPAR}

#printL "`cat ${OUTDIR}/mapped_reads_annotation_summary.txt`"
printL "`cat ${annotSummary}`"

if [ "${isBAM}" = 1 ]; then
  # insilico cut stats
  module="module1"
  step="insilico_cut"
  outRscript="${OUTDIR}/${module}_${step}.r"
  fprefix=`basename "${INFILE}"`
  cat ${SPARPATH}/scripts/R/${module}/${step}/*.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'\""; print wdir; fprefix="fprefix=\"'${fprefix}'\""; print fprefix; maxReadLength="maxReadLength='${maxReadLength}'"; print maxReadLength;}{print}' - > ${outRscript}
  ${RSCRIPT} ${outRscript}
fi

if [ "${isFASTQ}" = 1 ]; then
  # mapping stats
  module="module1"
  step="mapping"
  outRscript="${OUTDIR}/${module}_${step}.r"
  submoduleDir=${SPARPATH}/scripts/R/${module}/${step}
  cat ${submoduleDir}/*.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'\""; print wdir;}{print}' - > ${outRscript}
  ${RSCRIPT} ${outRscript}
  runScript "make_submodule_html.sh ${submoduleDir}" > ${OUTDIR}/${module}_${step}.html  

  # insilico cut stats
  module="module1"
  step="insilico_cut"
  outRscript="${OUTDIR}/${module}_${step}.r"
  fprefix=`basename "${OUTBAM}"`
  submoduleDir=${SPARPATH}/scripts/R/${module}/${step}
  cat ${submoduleDir}/*.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'\""; print wdir; fprefix="fprefix=\"'${fprefix}'\""; print fprefix; maxReadLength="maxReadLength='${maxReadLength}'"; print maxReadLength;}{print}' - > ${outRscript}
  ${RSCRIPT} ${outRscript}
  runScript "make_submodule_html.sh ${submoduleDir}" > ${OUTDIR}/${module}_${step}.html  
fi


