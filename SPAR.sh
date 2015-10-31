

source `dirname $0`/config.sh

set -e

INFASTQ=$1

if [ $# -lt 1 ]
then
  echo "USAGE: `basename $0` reads.fastq <SPAR_output_directory>"
  exit 1
fi


if [ ! -f "${INFASTQ}" ]; then
  echo -e "*****ERROR: FASTQ file\n${INFASTQ}\nnot found!"
  exit 1
fi


TRIMFASTQ=${INFASTQ}


EXPNAME=`basename ${TRIMFASTQ}`
EXPNAME=${EXPNAME%_*}

SPARDIR=${2:-${SPARDIR}}
SPARPATH=`dirname $0`

OUTDIR=${SPARDIR}/${EXPNAME}_m${maxMismatchCnt}_map${maxMapCnt}
mkdir -p ${OUTDIR}

LOGSPAR=${OUTDIR}/SPAR.log
>${LOGSPAR}

OUTBAM=${OUTDIR}/Aligned.out.filtered.hardClipped.sorted.bam

#echo "EXPNAME=${EXPNAME}"
#echo "${OUTBAM}"

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

printT "SPAR run started"

runScript "run_star_smrna2.sh ${TRIMFASTQ} ${maxMismatchCnt} ${maxMapCnt} ${OUTDIR}"

printT "Converting BAM to bedGraph"
runScript "bam_to_bedgraph.sh ${OUTBAM}"

printT "Converting bedGraph to bigWig"
runScript "bedgraph_to_bigwig.sh ${OUTBAM}.pos.bedgraph"
runScript "bedgraph_to_bigwig.sh ${OUTBAM}.neg.bedgraph"

printT "Segmenting [positive strand]"
runScript "segment_bedgraph_entropy.sh ${OUTBAM}.pos.bedgraph pos"
printT "Segmenting [negative strand]"
runScript "segment_bedgraph_entropy.sh ${OUTBAM}.neg.bedgraph neg"

printT "Annotating [positive strand]"
runScript "annotate_segm2.sh ${OUTBAM}.pos.bedgraph.segm"
printT "Annotating [negative strand]"
runScript "annotate_segm2.sh ${OUTBAM}.neg.bedgraph.segm"


printT "DONE."


#chr,a,b,peakID,exprVal,strand,entropy5p,normEntropy5p,maxEntropy5p, entropy3p, normEntropy3p, maxEntropy3p, maxProb5p, maxProb3p;

# add headers
awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction""annotChr", "annotChrStart", "annotChrEnd", "annotID", "annotRNAclass","annotStrand","overlap"; exit}' |  cat - ${OUTBAM}.*.bedgraph.segm.annot.final > ${OUTBAM}.annot.final
awk 'BEGIN{OFS="\t"; print "#chr", "chrStart", "chrEnd", "peakID", "expressionValue", "strand", "entropy5p", "normalizedEntropy5p", "maxEntropy5p","entropy3p","normalizedEntropy3p","maxEntropy3p","same5pFraction","same3pFraction"; exit}' |  cat - ${OUTBAM}.*.bedgraph.segm.unannotated.bed > ${OUTBAM}.unannot
#cat ${OUTBAM}.*.bedgraph.segm.annot.final > ${OUTBAM}.annot.final
#cat ${OUTBAM}.*.bedgraph.segm.unannotated.bed > ${OUTBAM}.unannot

awk 'BEGIN{OFS="\t";}{if ($0~/^#/) {next}; if (NR==FNR) {exprVal=$5; rnaClass=$19; exprPerClass[rnaClass]+=exprVal;classCnt[rnaClass]+=1;totalExprAnnot+=exprVal; totalAnnotPeakCnt+=1}else{exprVal=$5; totalExprUnannot+=exprVal;totalUnannotPeakCnt+=1;}}END{totalPeakCnt=totalAnnotPeakCnt+totalUnannotPeakCnt; totalExpr=totalExprAnnot+totalExprUnannot; for (rnaClass in exprPerClass) print rnaClass, classCnt[rnaClass], exprPerClass[rnaClass], exprPerClass[rnaClass]/totalExpr; print "Annotated", totalAnnotPeakCnt, totalExprAnnot, totalExprAnnot/totalExpr; print "Unannotated",totalUnannotPeakCnt,totalExprUnannot,totalExprUnannot/totalExpr}' ${OUTBAM}.annot.final ${OUTBAM}.unannot | sort -k1,1 | awk 'BEGIN{OFS="\t"; print "#RNA class","Peaks","Reads","Fraction of reads"}{print}' > ${OUTDIR}/mapped_reads_annotation_summary.txt


printL "\n===Output==="
printL "Output directory: ${OUTDIR}"
printL "\nMapping output:"
printL "${OUTBAM}"

printL "Annotation output:"
#ls ${OUTBAM}.*.bedgraph.segm.annot.final | tee -a ${LOGSPAR}
ls ${OUTBAM}.annot.final | tee -a ${LOGSPAR}

printL "Un-annotated output:" 
#ls ${OUTBAM}.*.bedgraph.segm.unannotated.bed | tee -a ${LOGSPAR}
ls ${OUTBAM}.unannot | tee -a ${LOGSPAR}


printL "\nAnnotation summary:"
ls ${OUTDIR}/mapped_reads_annotation_summary.txt | tee -a ${LOGSPAR}

printL "\nMapping stats:"
ls ${OUTDIR}/MAPSTAT.txt | tee -a ${LOGSPAR}

printL "\nMapped reads alignment summary:"
ls ${OUTDIR}/cigar.stat | tee -a ${LOGSPAR}


printL "\n\n===Run summary==="
printL "FASTQ: ${TRIMFASTQ}"
grep -e "Reads \[all\]" ${OUTDIR}/MAPSTAT.txt | awk 'BEGIN{FS="\t"}{printf "Mapped reads: %d [%.4f%%]\n", $2, $3}' | tee -a ${LOGSPAR}

numAnnot=$(cat ${OUTBAM}.*.bedgraph.segm.annot.final | wc -l)
printL "Annotated loci count: ${numAnnot}"
#printL "Annotated loci by RNA class:"
#cat ${OUTBAM}.*.bedgraph.segm.annot.final | cut -f 19 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$1}' | sort -k2,2nr | tee -a ${LOGSPAR}

numUnannot=$(cat ${OUTBAM}.*.bedgraph.segm.unannotated.bed | wc -l)
printL "Un-annotated loci count: ${numUnannot}"


printL "\nAnnotation summary:"
ls ${OUTDIR}/mapped_reads_annotation_summary.txt | tee -a ${LOGSPAR}

printL "`cat ${OUTDIR}/mapped_reads_annotation_summary.txt`"
