plot_ggplot2_TRIMMING_STATS<- function(datafile="SRR527595.fastq.trimming.stepC.stats",wdir=".") {

# Plots for Module 1
# Trimming and Mapping statistics 

# Module_1_Figure_1 (Figure 1.1) 
# description of figure: No of reads of that are 3' and 5' trimmed, 3' trimmed only, 5' trimmed only or have no adapters at all
# input: sampleID.fastq.trimming.stepC.stats
# output: TRIMMING_STATS.png / TRIMMING_STATS.pdf
# comment: TRIMMING stats - number of reads that have both 3' and 5' adapters, 3' adapters only, 5' adapters only, or no adapters 

BASENAME="TRIMMING_STATS"
PLOTTITLE="Trimming statistics"
XTITLE="Category"
YTITLE="Number of reads"

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))

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
DN = D[,c(1,3)]
DN$V1=c("3' and 5' trimmed","3' trimmed","5' trimmed","No adapters")
DN$V3 = as.numeric(DN$V3)
DN$V2 = DN$V3/sum(DN$V3) * 100
DN_f = ddply(DN, .(V1), transform, pos = (cumsum(V3) + 0.05*max(DN[,2])))
DN_f$label = paste0(sprintf("%.1f", DN_f$V2), "%")
colnames(DN_f) = c(XTITLE,YTITLE,"Percentage","pos","label")

options(scipen=10000)
chart = ggplot(DN_f, aes(x = DN_f[,XTITLE], y = DN_f[,YTITLE], fill = Percentage)) +
   geom_bar(stat = "identity", width = .7, fill="blue") +
   geom_text(aes(y = pos, label = label), size = 5,color="red")+
   xlab(XTITLE)+ylab(YTITLE)+ggtitle(PLOTTITLE)+theme_grey(base_size = 12)+
   theme(legend.position="none",axis.text = element_text(size = 12),axis.title = element_text(size =14),plot.title = element_text(size = 16))

plot(chart)
ggsave(pdffile,dpi=600)
ggsave(pngfile,dpi=600)

}

plot_ggplot2_TRIMMING_STATS(wdir=wdir)
