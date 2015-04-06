#Code for plotting outputs from Metagenomic analysis
#NOTE: Please only use this code after running the EENB_Analysis file to get csv output!
readFiles <- function(dir) {
  setwd(dir)
  files <- (Sys.glob("*.csv"))
  listOfFiles <- lapply(files, function(x) read.table(x, sep=",", header=TRUE))
  return(listOfFiles)
}
analyses<-readFiles("~/Metagenomic Data/Full Outputs")
#Calculate effect sizes between watered/drought groups and rank them
#then plot effect sizes by mean groups for visual of significant differences
for (i in 1:length(analyses)){
  analyses[[i]]$effect<-analyses[[i]]$mean_group1-analyses[[i]]$mean_group2
  analyses[[i]]$rank<-rank(analyses[[i]]$effect, ties.method="average")
}

#Create some plots to summarize the data
#Can remove the geom_text option to remove rankings
library(ggplot2)
source("multiplot_code.R")
data.names<-c("Dataset ECIDs.cts_75", "Dataset ECIDs.cts_75 (method 2)",
              "Dataset GOProcIDs.cts_75", "Dataset GOProcIDs.cts_75 (method 2)",
              "Dataset GOFuncIDs.cts_75","Dataset GOFuncIDs.cts_75 (method 2)",
              "Dataset PFamIDs.cts_75","Dataset PFamIDs.cts_75 (method 2)", 
              "Dataset TIGRFamIDs.cts_75","Dataset TIGRFamIDs.cts_75 (method 2)","Dataset Phyllo_Genus.cts_75","Dataset Phyllo_Genus.cts_75 (method 2)")
outnames<-paste(paste("Summary Plots", data.names, sep="-"),".pdf",sep="")
#Outnames not working in pdf file, so for now have to enter file names manually...
#Will look into this for future use
#For now, just change the index of i and the name to the outputs manually
i<-1
pdf("Summary Plots-Dataset Phyllo_Genus.cts_75 (method 2).pdf")
ggplot(analyses[[i]], aes(x=analyses[[i]]$mean_group1,y=analyses[[i]]$effect, colour=analyses[[i]]$effect))+geom_point(size=3.5)+
  labs(x="Mean: Watered Group", y="Effect Size")+scale_colour_gradient("Effect",low="green", high="blue")+
  ggtitle(expression(atop("Mean of Watered Group by Effect Size")))+scale_x_reverse()+
  theme(plot.title = element_text(size=20, face="bold", colour="black", vjust=-1), axis.title=element_text(size=14))

ggplot(analyses[[i]], aes(x=analyses[[i]]$mean_group2,y=analyses[[i]]$effect, colour=analyses[[i]]$effect))+geom_point(size=3.5)+
  labs(x="Mean: Drought Group", y="Effect Size")+scale_colour_gradient("Effect",low="red", high="yellow")+
  ggtitle(expression(atop("Mean of Drought Group by Effect Size")))+scale_x_reverse()+
  theme(plot.title = element_text(size=20, face="bold", colour="black", vjust=-1), axis.title=element_text(size=14))
  
ggplot(analyses[[i]], aes(x=analyses[[i]]$effect,y=-log10(analyses[[i]]$p.adj), colour=-log10(analyses[[i]]$p.adj)))+geom_point(size=3.5)+
  labs(x="Effect Size", y="-log10(p-value)")+scale_colour_gradientn("Transformed\np-value",colours=rainbow(5))+
  ggtitle(expression(atop("Effect Size by -log10 Adjusted p-value")))+scale_x_reverse()+
  theme(plot.title = element_text(size=20, face="bold", colour="black", vjust=-1), axis.title=element_text(size=14))
dev.off()

#Higher res SVG plots
svg("Watered by Effect-Dataset ECIDs.cts_75.svg")
ggplot(analyses[[i]], aes(x=analyses[[i]]$mean_group1,y=analyses[[i]]$effect, colour=analyses[[i]]$effect))+geom_point(size=3.5)+
  labs(x="Mean: Watered Group", y="Effect Size")+scale_colour_gradient("Effect",low="green", high="blue")+
  ggtitle(expression(atop("Mean of Watered Group by Effect Size")))+scale_x_reverse()+
  theme(plot.title = element_text(size=20, face="bold", colour="black", vjust=-1), axis.title=element_text(size=14))
dev.off()

svg("Drought by Effect-Dataset ECIDs.cts_75.svg")
ggplot(analyses[[i]], aes(x=analyses[[i]]$mean_group2,y=analyses[[i]]$effect, colour=analyses[[i]]$effect))+geom_point(size=3.5)+
  labs(x="Mean: Drought Group", y="Effect Size")+scale_colour_gradient("Effect",low="red", high="yellow")+
  ggtitle(expression(atop("Mean of Drought Group by Effect Size")))+scale_x_reverse()+
  theme(plot.title = element_text(size=20, face="bold", colour="black", vjust=-1), axis.title=element_text(size=14))
dev.off()

svg("Effect by Pvalue-Dataset ECIDs.cts_75.svg")
ggplot(analyses[[i]], aes(x=analyses[[i]]$effect,y=-log10(analyses[[i]]$p.adj), colour=-log10(analyses[[i]]$p.adj)))+geom_point(size=3.5)+
  labs(x="Effect Size", y="-log10(p-value)")+scale_colour_gradientn("Transformed\np-value",colours=rainbow(5))+
  ggtitle(expression(atop("Effect Size by -log10 Adjusted p-value")))+scale_x_reverse()+
  theme(plot.title = element_text(size=20, face="bold", colour="black", vjust=-1), axis.title=element_text(size=14))
dev.off()