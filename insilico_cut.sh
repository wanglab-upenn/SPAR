INBAM=$1
maxReadLength=$2
SAMTOOLS=samtools
CUTBAM=${INBAM%.*}.max${maxReadLength}.bam
OTHERBAM=${INBAM%.*}.greaterthan${maxReadLength}.bam
${SAMTOOLS} view -h ${INBAM} | \
awk 'BEGIN{FS="\t"; OFS="\t"; maxLength='${maxReadLength}'+0; otherFile="'${otherFile}'"; samtoolsCut="'${SAMTOOLS}' view -bS - > '${CUTBAM}'"; samtoolsOther="'${SAMTOOLS}' view -bS - > '${OTHERBAM}'";}
     {
        if ($0~/^@/) { print | samtoolsCut; print $0 | samtoolsOther; next } # print header

        cigar=$6;
        if (cigar~/^[0-9]+M$/)
        {
            
            n=split(cigar,a,"M");
            readLength=a[1];
            if (readLength<=maxLength)
            {
              # output genomic alignments with length \leq maxLength
              print | samtoolsCut;
            }
            else
            {
              # output alignments with length > maxLength
              print | samtoolsOther;
            }
        }
        else
        {
          # print other alignments
          print | samtoolsOther # > otherFile 
        }
     }' #| \
#${SAMTOOLS} view -bS - > ${CUTBAM} 

 
