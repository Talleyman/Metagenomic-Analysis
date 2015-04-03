#Code for plotting outputs from Metagenomic analysis
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
library(ggplot2)