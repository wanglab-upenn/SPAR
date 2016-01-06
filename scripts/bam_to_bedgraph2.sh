set -e
source `dirname $0`/../config.sh
INBAM=$1

BAM=`basename ${INBAM}`
BAMDIR=`dirname ${INBAM}`


# directory where output will be saved
# directory is either provided as input argument or is the same as the input BAM directory
OUTDIR=${2:-${BAMDIR}}
mkdir -p ${OUTDIR}

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
    BEGIN{OFS="\t";chrLen=0;verbose=0;}
     {
       if (NR==1) chr_prev = $3;
       chr = $3
       if (chr == chr_prev)
       {
         flag = $2 # SAM flag
         strand = "+";
         if (and(flag,16)>0) strand="-";
         rstart = $4+0 # read start
         rend = $4+length($10)-1+0 # read end
         nHits = substr($12,6);
         w = 1 / nHits;
         if (chrLen<rend) chrLen=rend;
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
     }' 

# need to sort because of unnatural order of samtools view 
# sort by chr, chrStart, chrEnd
sort -k1,1 -k2,2n -k3,3n ${OUTDIR}/${BAM}.pos.bedgraph -o ${OUTDIR}/${BAM}.pos.bedgraph
sort -k1,1 -k2,2n -k3,3n ${OUTDIR}/${BAM}.neg.bedgraph -o ${OUTDIR}/${BAM}.neg.bedgraph

printT "BEDGRAPH finished successfuly."
