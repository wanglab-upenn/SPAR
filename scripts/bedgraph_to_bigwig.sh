set -e
bgfile=$1

source `dirname $0`/../config.sh
outbigwig=${bgfile/bedgraph/bigWig}
${BGTOBIGWIG} ${bgfile} ${chromInfo} ${outbigwig}
