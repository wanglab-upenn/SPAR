set -e
bgfile=$1

source `dirname $0`/../config.sh
outbigwig=${bgfile/bedgraph/bigWig}
chromInfo=`dirname $0`/../annot/chromInfo.txt
${BGTOBIGWIG} ${bgfile} ${chromInfo} ${outbigwig}
