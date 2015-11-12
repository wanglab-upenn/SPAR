#SPAR config file

HOMEDIR=/home/pkuksa #${HOME}

#absolute path to the bin directory
BINDIR=${HOMEDIR}/bin

#absolute path to the SPAR home directory
#SPARPATH=${BINDIR}/SPAR_github/SPAR

# (default) absolute path to the SPAR output directory
SPARDIR=${HOMEDIR}/SPAR_out

#absolute path to pre-installed STAR, samtools, AWK, etc

# MAPPING
STAR=${BINDIR}/STAR-STAR_2.4.0k/bin/Linux_x86_64/STAR # STAR

SAMTOOLS=${BINDIR}/samtools-1.2/samtools # SAMTOOLS
#BEDTOOLS=${BINDIR}/bedtools2/bin/bedtools # BEDTOOLS
BEDTOOLS=/usr/local/bin/bedtools # BEDTOOLS

GAWK=${BINDIR}/gawk-4.1.0/gawk

# UCSC tools
BGTOBIGWIG=${BINDIR}/bedGraphToBigWig
BEDTOBIGBED=${BINDIR}/bedToBigBed

# Adapter trimming
CUTADAPT=${BINDIR}/cutadapt-1.8.1/bin/cutadapt

#absolute path to the STAR genome index
genomeDir=${HOMEDIR}/datasets/hg19/star_2.4.0k/  # STAR genome index

#hg 19 chromosome information file
#chromInfo=${SPARPATH}/annot/chromInfo.txt

#mapping parameters for STAR
maxMismatchCnt=0 # maximum number of genomic mismatches
maxMapCnt=100 # maximum number of places a read can map to
minMappedLength=14 # minimum *mapped* length


function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1"
}
