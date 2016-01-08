plot_ggplot2_CLIPPING_STATS<- function(datafile="MAPSTAT.txt",wdir=".") {

# Plots for Module 1
# Trimming and Mapping statistics 

# Module_1_Figure_4 (Figure 1.4) 
# description of figure: Total number of mapped reads, reads with 5'-clipping, reads with no soft clipping and reads with soft clipping
# input: MAPSTAT.txt.forR
# output: CLIPPING_STATS.png / CLIPPING_STATS.pdf
# comment: Clipping stats - number of reads that are with 5'clipping, soft-clipping and without soft-clipping 

datafile=paste(wdir, "/", datafile, sep="")
BASENAME="CLIPPING_STATS"
PLOTTITLE="Clipping statistics on reads <44nt"
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
DN = D[c(4,11,12),-3]
DN[4,] = c("Five", DN$V2[1]-DN$V2[2]-DN$V2[3])
DN$V1=c("Total","No soft clipping","With soft clipping","5' clipped")
DN$V2 = as.numeric(DN$V2)
DN=DN[-1,]
DN$V3 = DN$V2/sum(DN$V2) * 100

DN_f = ddply(DN, .(V1), transform, pos = (cumsum(V2) - 0.5 * V2))
DN_f$label = paste0(sprintf("%.1f", DN_f$V3), "%")
colnames(DN_f) = c(XTITLE,YTITLE,"Percentage","pos","label")

options(scipen=10000)

chart = ggplot(DN_f, aes(x = DN_f[,XTITLE], y = DN_f[,YTITLE], fill = Percentage)) +
   geom_bar(stat = "identity", width = .7, fill="blue") +
   geom_text(aes(y = pos, label = label), size = 5,color="pale green")+
   xlab(XTITLE)+ylab(YTITLE)+ggtitle(PLOTTITLE)+theme_grey(base_size = 12)+
   theme(legend.position="none",axis.text = element_text(size = 12),axis.title = element_text(size =14),plot.title = element_text(size = 16))

plot(chart)
ggsave(pdffile,dpi=600)
ggsave(pngfile,dpi=600)

}

plot_ggplot2_CLIPPING_STATS(wdir=wdir)
