set -e
source `dirname $0`/../config.sh
INBAM=$1

BAM=`basename ${INBAM}`
BAMDIR=`dirname ${INBAM}`


# directory where output will be saved
# directory is either provided as input argument or is the same as the input BAM directory
OUTDIR=${2:-${BAMDIR}}
mkdir -p ${OUTDIR}

maxReadLength=${3:-44}


${SAMTOOLS} view ${INBAM} | \
${GAWK} 'function cmp_num_idx(i1, v1, i2, v2)
      {
       # numerical index comparison, ascending order
       return (i1 - i2)
      }
      function process_chr_coverage( c, n, chr, outfile)
      {
         bs = 1; # block start
         cnt = 0;
         f = 0; # read coverage
         f_prev = 0;
         PROCINFO["sorted_in"]="@ind_num_asc"
         for (i in c)
         {
            ++cnt;
            if (verbose==1 && cnt % 100000 == 0) printf "Processed=%d\n",cnt;
            #if (cnt==1) {bs=i; f_prev = c[i];}
            f+=c[i] # coverage
            # print i
            # print i, c[i], f
            if (f != f_prev)
            {
             #if (f_prev > 0.0001)
             # output intervals will be in UCSC format [0-based, half-open) 
             if (f_prev >= 1)
               #print chr, bs-1, i-1, f_prev > outfile
               printf "%s\t%d\t%d\t%.4f\n", chr, bs-1, i-1, f_prev > outfile
             bs = i;
            }
            f_prev = f;
         }
      }
    BEGIN{OFS="\t";chrLen=0;verbose=0;maxReadLength='${maxReadLength}'+0;}
     {
       if (NR==1) chr_prev = $3;
       chr = $3
       if (chr == chr_prev)
       {
         flag = $2 # SAM flag
         strand = "+";
         if (and(flag,16)>0) strand="-";
         rstart = $4+0 # read start
         rlength=length($10) # read length
         rend = $4+rlength-1+0 # read end
         nHits = substr($12,6);
         w = 1 / nHits;
         cigar=$6;
         if (chrLen<rend) chrLen=rend;

         if (rlength > maxReadLength)
         {
           cxh=(cigar "\t" nHits);
           aboveMaxReadCnt+=w;
           aboveCigarCount[cigar]+=1;
           aboveCigarCountWeighted[cigar]+=w;
           aboveCigarXnumHits[cxh]++;
           if (nHits == 1)
             aboveCigarCountUniq[cigar]++;
           # insilico cut: skip reads with length above max
           next;
         }

         # compute CIGAR stats
         cxh=(cigar "\t" nHits);
         belowCigarCount[cigar]++;
         belowCigarCountWeighted[cigar]+=w; 
         belowCigarXnumHits[cxh]++;
         if (nHits == 1)
           belowCigarCountUniq[cigar]++;
         
         # number of ~reads
         belowMaxReadCnt+=w;

         # print $0
         # print strand, rstart, rend, nHits, w
         if (strand=="+")
         {
             rcov_plus[rstart]+=w;
             rcov_plus[rend+1]-=w;
         }
         else
         {
             rcov_minus[rstart]+=w;
             rcov_minus[rend+1]-=w;
         } 
         if ( verbose==1 && (NR % 1000000) == 0 ) printf "Processed %d read alignments\n", NR;
       }
       else # process chromosome
       {
         n = chrLen; #chrLen[chr_prev]
         if (verbose==1) 
           printf "Processing %s [len=%d, Watson strand]\n", chr_prev, n
         outfile = "'${OUTDIR}'/'${BAM}'.pos.bedgraph"
 
         process_chr_coverage( rcov_plus, n, chr_prev, outfile )     
         split("",rcov_plus,":") # clear coverage array
         if (verbose==1)
           printf "Processing %s [len=%d, Crick strand]\n", chr_prev, n
         outfile = "'${OUTDIR}'/'${BAM}'.neg.bedgraph"
         process_chr_coverage( rcov_minus, n, chr_prev, outfile )     
         split("",rcov_minus,":") # clear coverage array
         chrLen = 0;
       }
       chr_prev = chr;        
     }
     END{
         # process last chromosome
         n = chrLen; #chrLen[chr_prev]
         if (verbose==1)
           printf "Processing %s [len=%d, Watson strand]\n", chr_prev, n
         outfile = "'${OUTDIR}'/'${BAM}'.pos.bedgraph"
         process_chr_coverage( rcov_plus, n, chr_prev, outfile )     
         split("",rcov_plus,":") # clear coverage array
         if (verbose==1)
           printf "Processing %s [len=%d, Crick strand]\n", chr_prev, n
         outfile = "'${OUTDIR}'/'${BAM}'.neg.bedgraph"
         process_chr_coverage( rcov_minus, n, chr_prev, outfile )     
         split("",rcov_minus,":") # clear coverage array
         chrLen = 0;
 
         printf "max%d:\t%d\ngreaterthan%d:\t%d\n", maxReadLength, belowMaxReadCnt, maxReadLength, aboveMaxReadCnt > "'${OUTDIR}'/'${BAM}'.insilico_cut.stats"

         # output CIGAR stats for reads above the max length 
         cigarstatsFile=("'${OUTDIR}'/'${BAM}'.greaterthan" maxReadLength ".cigar.stats")
         for (c in aboveCigarCount)
         {
               printf "%s\t%d\t%9.0f\t%d\n", c, aboveCigarCount[c], aboveCigarCountWeighted[c], aboveCigarCountUniq[c] > cigarstatsFile

         }
         cigarXnumHitsFile=("'${OUTDIR}'/'${BAM}'.greaterthan" maxReadLength ".cigar.nhits.stats")
         for (c in aboveCigarXnumHits)
         {
            printf "%s\t%d\n", c, aboveCigarXnumHits[c] > cigarXnumHitsFile 
         }
 
         # output CIGAR stats for reads at or below max length 
         cigarstatsFile=("'${OUTDIR}'/'${BAM}'.max" maxReadLength ".cigar.stats")
         for (c in belowCigarCount)
         {
               printf "%s\t%d\t%9.0f\t%d\n", c, belowCigarCount[c], belowCigarCountWeighted[c], belowCigarCountUniq[c] > cigarstatsFile

         }
         cigarXnumHitsFile=("'${OUTDIR}'/'${BAM}'.max" maxReadLength ".cigar.nhits.stats")
         for (c in belowCigarXnumHits)
         {
            printf "%s\t%d\n", c, belowCigarXnumHits[c] > cigarXnumHitsFile 
         }
 
         
     }' 

# need to sort because of unnatural order of samtools view 
# sort by chr, chrStart, chrEnd
sort -k1,1 -k2,2n -k3,3n ${OUTDIR}/${BAM}.pos.bedgraph -o ${OUTDIR}/${BAM}.pos.bedgraph
sort -k1,1 -k2,2n -k3,3n ${OUTDIR}/${BAM}.neg.bedgraph -o ${OUTDIR}/${BAM}.neg.bedgraph

printT "BEDGRAPH finished successfuly."
