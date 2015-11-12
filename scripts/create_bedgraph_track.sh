set -e
source `dirname $0`/../config.sh

if [ $# -lt 1 ]
then
  echo "USAGE: `basename $0` .bedgraph"
  exit 1
fi

chromInfo=`dirname $0`/../annot/chromInfo.txt
INBG=$1 # input bedgraph
OUTBIGWIG="${INBG%.*}.bigWig" # output bedgraph
${BGTOBIGWIG} ${INBG} ${chromInfo} ${OUTBIGWIG}

#bedGraphToBigWig
