plot_ggplot2_FASTQ_READ_LENGTH <- function(datafile="SRR527595_A1_B2_C.fastq.read_length",wdir=".") {

# Plots for Module 1
# Trimming and Mapping statistics 

# Module_1_Figure_2 (Figure 1.2) 
# description of figure: fastq read length figure 
# input: sampleID.fastq.read_length, e.g. SRR527595_A1_B2_C.fastq.read_length
# output: FASTQ_READ_LENGTH.png / FASTQ_READ_LENGTH.pdf
# comment: length distribution of trimmed reads 

BASENAME="FASTQ_READ_LENGTH"
PLOTTITLE="Length distribution of FASTQ reads"
XTITLE="FASTQ read length"
YTITLE="Count"

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(png))
#suppressPackageStartupMessages(library(Cairo))

#args<-commandArgs(TRUE)
#datafile=args[1] # input file
#wdir=args[2] # output / working directory
#if (length(args)<1) { stop("ERROR: No input! USAGE: script inputfile <output-dir>")}
#if (length(args)<2) { wdir="." } 
#use current dir if no 
#working dir has been specified

# output image file
pngfile= paste(wdir, "/", paste(BASENAME,".png",sep=""), sep="")
pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")
# output html file
htmlfile=paste(wdir, "/", paste(BASENAME,".html",sep=""), sep="")

D = read.table(datafile,header=F,sep='\t')
i1 = grep("D",D$V1)
u1 = D[i1,]
'%!in%' <- function(x,y)!('%in%'(x,y))
x = D[D$V1%!in%u1$V1,]
i2 = grep("I",x$V1)
u2 = x[i2,]
f = x[x$V1%!in%u2$V1,]
f = f[,-3]
xf <- lapply(f, function(x) {gsub("M", "", x)})
xf <- matrix(unlist(xf), ncol = 2)
colnames(xf)=c("V1","A1B2")
rownames(xf)=xf[,1]
xf = apply(xf, 2, as.numeric)
xf = data.frame(xf)
colnames(xf) = c(XTITLE, YTITLE)

options(scipen=10000)
chart = ggplot(xf, aes(x=xf[,XTITLE], y=xf[,YTITLE]))+  
stat_summary(fun.y=mean,geom="bar", aes(width=0.2))+
geom_bar(stat = "identity", fill="blue") +
ggtitle(PLOTTITLE) + xlab(XTITLE)+ylab(YTITLE) +
theme(legend.position="none",axis.text = element_text(size = 12),axis.title = element_text(size =14),plot.title = element_text(size = 16))

plot(chart)
ggsave(pdffile,dpi=600)
ggsave(pngfile,dpi=600, type="cairo-png")

}

plot_ggplot2_FASTQ_READ_LENGTH(wdir=wdir)

