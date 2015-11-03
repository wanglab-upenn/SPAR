#!/bin/bash

set -e

#source `dirname $0`/../config.sh
CUTADAPT=${HOME}/bin/cutadapt-1.9/bin/cutadapt

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


${CUTADAPT} ${ADAPTERoption} -e $MAX_ADAPTER_ERROR -o $TRIMMEDREADS \
--untrimmed-output $UNTRIMMEDREADS --too-short-output $TOOSHORTREADS \
-O $MIN_ADAPTER_MATCH -m $MIN_TRIMMED_READ_LEN -n 2 \
--length-tag "length=" $RAWREADS


#${CUTADAPT} -a ${ADAPTER3p} -g ${ADAPTER5p} -e $MAX_ADAPTER_ERROR -o $TRIMMEDREADS \
#--untrimmed-output $UNTRIMMEDREADS --too-short-output $TOOSHORTREADS \
#-O $MIN_ADAPTER_MATCH -m $MIN_TRIMMED_READ_LEN -n 2 \
#--length-tag "length=" $RAWREADS
