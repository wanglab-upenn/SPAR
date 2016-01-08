#plot_ggplot2_MAPPED_READ_LENGTH_INSILICO_COMP_UNIQvsMULTI <- function(datafile1="Aligned.out.filtered.hardClipped.sorted.bam.max44.cigar.stats",datafile2="Aligned.out.filtered.hardClipped.sorted.bam.greaterthan44.cigar.stats",wdir=".",maxReadLength=44) {
plot_ggplot2_MAPPED_READ_LENGTH_INSILICO_COMP_UNIQvsMULTI <- function(fprefix="Aligned.out.filtered.hardClipped.sorted.bam",wdir=".",maxReadLength=44) {

datafile1=paste(wdir, "/", fprefix, ".max", maxReadLength, ".cigar.stats",sep="")
datafile2=paste(wdir, "/", fprefix, ".greaterthan", maxReadLength, ".cigar.stats",sep="")

# Plots for Module 1
# Trimming and Mapping statistics 

# Module_1_Figure_5 (Figure 1.5) 
# description of figure: Mapped read length figure - characteristics of reads that are <=44nt, >44nt, uniquely and multi-mapped. 
# input: Aligned.out.filtered.hardClipped.sorted.max44.bam.cigar.stats, Aligned.out.filtered.hardClipped.sorted.greaterthan44.bam.cigar.stats
# output: MAPPED_READ_LENGTH_COMPARISON.png / MAPPED_READ_LENGTH_COMPARISON.pdf
# comment: length distribution of mapped reads (of different length, and whether they are uniquely- or multi-mapped) 

BASENAME="MAPPED_READ_LENGTH_INSILICO_COMP_UNIQvsMULTI"
PLOTTITLE1="Frequency of each read (<=44nt) being mapped"
PLOTTITLE2="Frequency of each read (>44nt) being mapped"
XTITLE="Mapped read length"
YTITLE="Count"

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(gridExtra))

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
f = f[,-2]
xf <- lapply(f, function(x) {gsub("M", "", x)})
xf <- matrix(unlist(xf), ncol = 3)
rownames(xf)=xf[,1]
xf = apply(xf, 2, as.numeric)
xf = data.frame(xf)
colnames(xf) = c("X1", "Multi-mapped","Uniquely-mapped")
xf = melt(xf,id.vars='X1')
colnames(xf) = c(XTITLE, "variable","value")
xf[,2] <- factor(xf[,2], levels =  c("Uniquely-mapped","Multi-mapped"))
xf = ddply(xf, c('variable'))
colnames(xf) = c(XTITLE, "variable",YTITLE)
D_max_xf = xf

D_greater = read.table(datafile2)
i1 = grep("D",D_greater$V1);u1 = D_greater[i1,]
'%!in%' <- function(x,y)!('%in%'(x,y))
x = D_greater[D_greater$V1%!in%u1$V1,]
i2 = grep("I",x$V1);u2 = x[i2,]
f = x[x$V1%!in%u2$V1,]
f = f[,-2]
xf <- lapply(f, function(x) {gsub("M", "", x)})
xf <- matrix(unlist(xf), ncol = 3)
rownames(xf)=xf[,1]
xf = apply(xf, 2, as.numeric)
xf = data.frame(xf)
colnames(xf) = c("X1", "Multi-mapped","Uniquely-mapped")
xf = melt(xf,id.vars='X1')
colnames(xf) = c(XTITLE, "variable","value")
xf[,2] <- factor(xf[,2], levels =  c("Uniquely-mapped","Multi-mapped"))
xf = ddply(xf, c('variable'))
colnames(xf) = c(XTITLE, "variable",YTITLE)
D_greater_xf = xf

options(scipen=10000)
chart1 = ggplot(D_max_xf,aes(x=D_max_xf[,XTITLE],y=D_max_xf[,YTITLE],fill=D_max_xf[,2]))+ 
geom_bar(stat = "identity")+ 
scale_fill_manual(values = c("red","blue"),guide = guide_legend(title = "")) +
# scale_y_continuous(limits = c(min(D_max_xf[,XTITLE]),max(D_greater_xf[,XTITLE]))) + 
ggtitle(PLOTTITLE1) + xlab(XTITLE)+ylab(YTITLE) +
theme(legend.position="right",legend.text=element_text(size=10),axis.text = element_text(size = 10),axis.title = element_text(size =10),plot.title = element_text(size = 12))

chart2 = ggplot(D_greater_xf,aes(x=D_greater_xf[,XTITLE],y=D_greater_xf[,YTITLE],fill=D_greater_xf[,2]))+ 
geom_bar(stat = "identity")+ 
scale_fill_manual(values = c("red","blue"),guide = guide_legend(title = "")) +
ggtitle(PLOTTITLE2) + xlab(XTITLE)+ylab(YTITLE) +
theme(legend.position="right",legend.text=element_text(size=10),axis.text = element_text(size = 10),axis.title = element_text(size =10),plot.title = element_text(size = 12))


#png(file = pngfile,width=6.38, height=6.36, units="in", res=600)
#print(plot(grid.arrange(chart1,chart2,nrow=2)))
#dev.off() 

#pdf(file = pdffile)
#print(plot(grid.arrange(chart1,chart2,nrow=2)))
#dev.off() 
chart3=grid.arrange(chart1,chart2,nrow=2)
ggsave(pdffile,chart3,dpi=600)
ggsave(pngfile,chart3,dpi=600)

# ggsave(pdffile,chart1,chart2,dpi=600)
# ggsave(pngfile,chart1,chart2,dpi=600)

}

plot_ggplot2_MAPPED_READ_LENGTH_INSILICO_COMP_UNIQvsMULTI(fprefix=fprefix,wdir=wdir,maxReadLength=maxReadLength)
