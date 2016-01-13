#plot_ggplot2_MAPPED_READ_LENGTH_FREQ_INSILICO_COMP <- function(datafile1="Aligned.out.filtered.hardClipped.sorted.max44.bam.cigar.nhits.stats",datafile2="Aligned.out.filtered.hardClipped.sorted.greaterthan44.bam.cigar.nhits.stats",wdir=".") {


plot_ggplot2_MAPPED_READ_LENGTH_FREQ_INSILICO_COMP <- function(fprefix="Aligned.out.filtered.hardClipped.sorted.bam",wdir=".",maxReadLength=44) {

datafile1=paste(wdir, "/", fprefix, ".max", maxReadLength, ".cigar.nhits.stats",sep="")
datafile2=paste(wdir, "/", fprefix, ".greaterthan", maxReadLength, ".cigar.nhits.stats",sep="")

# Plots for Module 1
# Trimming and Mapping statistics 

# Module_1_Figure_6 (Figure 1.6) 
# description of figure: Frequency of mapping distribution figure 
# input: Aligned.out.filtered.hardClipped.sorted.max44.bam.cigar.nhits.stats, Aligned.out.filtered.hardClipped.sorted.greaterthan44.bam.cigar.nhits.stats
# output: FREQ_MAPPED_DISTRIBUTION.png / FREQ_MAPPED_DISTRIBUTION.pdf
# comment: Frequency of a read being mappped

BASENAME="MAPPED_READ_LENGTH_FREQ_INSILICO_COMP"
PLOTTITLE1="Frequency of each read (<=44nt) being mapped"
PLOTTITLE2="Frequency of each read (>44nt) being mapped"
XTITLE="Number of reads"
YTITLE="Frequency"

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(gridExtra))
#suppressPackageStartupMessages(library(png))
#suppressPackageStartupMessages(library(Cairo))

#args<-commandArgs(TRUE)
#datafile1=args[1] # input file 1
#datafile2=args[2] # input file 2
#wdir=args[3] # output / working directory
#if (length(args)<1) { stop("ERROR: No input! USAGE: script inputfile <output-dir>")}
#if (length(args)<2) { stop("ERROR: No input! USAGE: script inputfile <output-dir>")}
#if (length(args)<3) { wdir="." } 
#use current dir if no 
#working dir has been specified

# output image file
pngfile= paste(wdir, "/", paste(BASENAME,".png",sep=""), sep="")
pdffile= paste(wdir, "/", paste(BASENAME,".pdf",sep=""), sep="")
# output html file
htmlfile=paste(wdir, "/", paste(BASENAME,".html",sep=""), sep="")

D_max = read.table(datafile1)
i1 = grep("D",D_max$V1);u1 = D_max[i1,]
'%!in%' <- function(x,y)!('%in%'(x,y))
x = D_max[D_max$V1%!in%u1$V1,]
i2 = grep("I",x$V1);u2 = x[i2,]
f = x[x$V1%!in%u2$V1,]
xf <- lapply(f, function(x) {gsub("M", "", x)})
xf <- matrix(unlist(xf), ncol = 3)
rownames(xf)=xf[,1]
xf = apply(xf, 2, as.numeric)
xf = data.frame(xf)
xf = ddply(xf,.(X2), summarize, sum=sum(X3))
xf[,3] = xf[,2]/xf[,1]
xf[,4] = cumsum(xf[,3])
colnames(xf) = c(YTITLE,"original","average", XTITLE) 
D_max_xf = xf

D_greater = read.table(datafile2)
i1 = grep("D",D_greater$V1);u1 = D_greater[i1,]
'%!in%' <- function(x,y)!('%in%'(x,y))
x = D_greater[D_greater$V1%!in%u1$V1,]
i2 = grep("I",x$V1);u2 = x[i2,]
f = x[x$V1%!in%u2$V1,]
xf <- lapply(f, function(x) {gsub("M", "", x)})
xf <- matrix(unlist(xf), ncol = 3)
rownames(xf)=xf[,1]
xf = apply(xf, 2, as.numeric)
xf = data.frame(xf)
xf = ddply(xf,.(X2), summarize, sum=sum(X3))
xf[,3] = xf[,2]/xf[,1]
xf[,4] = cumsum(xf[,3])
colnames(xf) = c(YTITLE,"original","average", XTITLE) 
D_greater_xf = xf

options(scipen=10000)
chart1 = ggplot(D_max_xf, aes(y=D_max_xf[,XTITLE], x=D_max_xf[,YTITLE]))+  
geom_point(size=3, color="blue") + 
scale_y_continuous(limits = c(0,max(D_max_xf[,XTITLE]))) + 
scale_x_continuous(limits = c(0,max(D_max_xf[,YTITLE]))) + 
ggtitle(PLOTTITLE1) + ylab(XTITLE)+xlab(YTITLE) +
theme(legend.position="right",legend.text=element_text(size=10),axis.text = element_text(size = 10),axis.title = element_text(size =10),plot.title = element_text(size = 12))

chart2 = ggplot(D_greater_xf, aes(y=D_greater_xf[,XTITLE], x=D_greater_xf[,YTITLE]))+  
geom_point(size=3, color="blue") + 
scale_y_continuous(limits = c(0,max(D_greater_xf[,XTITLE]))) +
scale_x_continuous(limits = c(0,max(D_max_xf[,YTITLE]))) +  
ggtitle(PLOTTITLE2) + ylab(XTITLE)+xlab(YTITLE) +
theme(legend.position="right",legend.text=element_text(size=10),axis.text = element_text(size = 10),axis.title = element_text(size =10),plot.title = element_text(size = 12))

chart3=grid.arrange(chart1,chart2,nrow=2)
ggsave(pdffile,chart3,dpi=600)
ggsave(pngfile,chart3,dpi=600,type="cairo-png")

#pdf(file = pdffile)
#print(plot(grid.arrange(chart1,chart2,nrow=2)))
#dev.off() 

#png(file = pngfile,width=6.38, height=6.36, units="in", res=600)
#print(plot(grid.arrange(chart1,chart2,nrow=2)))
#dev.off() 

# ggsave(pdffile,chart1,chart2,dpi=600)
# ggsave(pngfile,chart1,chart2,dpi=600)

}

plot_ggplot2_MAPPED_READ_LENGTH_FREQ_INSILICO_COMP(fprefix=fprefix,wdir=wdir,maxReadLength=maxReadLength) 
