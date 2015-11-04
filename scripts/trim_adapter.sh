#!/bin/bash

set -e

#source `dirname $0`/../config.sh
CUTADAPT=`which cutadapt 2> /dev/null | tail -n 1` || ${HOME}/bin/cutadapt-1.8.1/bin/cutadapt
echo ${CUTADAPT}
command -v ${CUTADAPT} > /dev/null 2>&1 || { echo >&2 "$0 requires cutadapt, but it is not found. Please cu"; exit 1; }

if [ $# -lt 1 ]
then
  echo "USAGE: `basename $0` reads.fastq <-a ADAPTERSEQ3p> <-g ADAPTERSEQ5p>"
  exit 1
fi



INFASTQ=$1
ADAPTERoption=$@
ADAPTERoption=${ADAPTERoption##${INFASTQ}}
#ADAPTER3p=$2
#ADAPTER5p=$3

DATASET=${INFASTQ%%.fastq}


# trimming parameters
MIN_TRIMMED_READ_LEN=14
MAX_ADAPTER_ERROR=0.06
MIN_ADAPTER_MATCH=6



TRIMMEDREADS=${DATASET}_trimmed.fastq
UNTRIMMEDREADS=${DATASET}_untrimmed.fastq
TOOSHORTREADS=${DATASET}_tooshort.fastq


${CUTADAPT} ${ADAPTERoption} -e ${MAX_ADAPTER_ERROR} -o ${TRIMMEDREADS} \
--untrimmed-output ${UNTRIMMEDREADS} --too-short-output ${TOOSHORTREADS} \
-O ${MIN_ADAPTER_MATCH} -m ${MIN_TRIMMED_READ_LEN} -n 2 \
--length-tag "length=" ${INFASTQ}


#${CUTADAPT} -a ${ADAPTER3p} -g ${ADAPTER5p} -e $MAX_ADAPTER_ERROR -o $TRIMMEDREADS \
#--untrimmed-output $UNTRIMMEDREADS --too-short-output $TOOSHORTREADS \
#-O $MIN_ADAPTER_MATCH -m $MIN_TRIMMED_READ_LEN -n 2 \
#--length-tag "length=" $RAWREADS
