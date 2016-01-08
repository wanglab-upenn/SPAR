plot_ggplot2_MAPPING_STATS <- function(datafile="MAPSTAT.txt",wdir=".") {

# Plots for Module 1
# Trimming and Mapping statistics 

# Module_1_Figure_3 (Figure 1.3) 
# description of figure: No of uniquely-, multi-, and un-mapped reads 
# input: MAPSTAT.txt.forR
# output: MAPPING_STATS.png / MAPPING_STATS.pdf
# comment: Mapping stats - number of uniquely-, multi-, and un-mapped reads 

datafile=paste(wdir, "/", datafile, sep="")

BASENAME="MAPPING_STATS"
#PLOTTITLE="Mapping statistics on reads <44nt"
PLOTTITLE="Read mapping statistics"
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
DN = D[1:3,-3]
DN[4,] = c("NA", DN$V2[1]-DN$V2[2]-DN$V2[3])
DN$V1=c("Total","Uniquely-mapped","Multi-mapped","Unmapped")
DN$V2 = as.numeric(DN$V2)
DN=DN[-1,]
DN$V3 = DN$V2/sum(DN$V2) * 100

DN_f = ddply(DN, .(V1), transform, pos = (cumsum(V2) - 0.5 * V2))
DN_f$label = paste0(sprintf("%.1f", DN_f$V3), "%")
colnames(DN_f) = c(XTITLE,YTITLE,"Percentage","pos","label")

options(scipen=10000)
print(DN_f)

print(class(DN_f[,YTITLE]))
#chart = ggplot(DN_f, aes(x = factor(DN_f[,XTITLE]) , y = factor(DN_f[,YTITLE]), fill = Percentage)) + 
chart = ggplot(DN_f, aes(x = DN_f[,XTITLE] , y = DN_f[,YTITLE], fill = Percentage)) + 
#scale_y_continuous(limits=c(0,100000000))+
#scale_y_discrete(breaks=seq(1,max(DN_f[,YTITLE])),0.1*max(DN_f[,YTITLE]))+
geom_bar(stat = "identity", width = .7, fill="blue") +
   geom_text(aes(y = DN_f[,"pos"], label = DN_f[,"label"]), size = 5,color="pale green")+
   #geom_text(aes(y = pos, label = label), size = 5,color="pale green")#+
   xlab(XTITLE)+ylab(YTITLE)+ggtitle(PLOTTITLE)+theme_grey(base_size = 12)+
   theme(legend.position="none",axis.text = element_text(size = 12),axis.title = element_text(size =14),plot.title = element_text(size = 16))

plot(chart)
ggsave(pdffile,dpi=600)
ggsave(pngfile,dpi=600)

}

plot_ggplot2_MAPPING_STATS(wdir=wdir)
