setwd("~/Metagenomic Data")
source("TwoStage_Package_Code.R")
#Read in data
ECID<-read.table("ECIDs.cts_75.Mapping.padded.main.cumulative.summary_table.tsv",sep="\t",quote="",header=TRUE)
ProcID<-read.table("GOProcIDs.cts_75.Mapping.padded.main.cumulative.summary_table.tsv",sep="\t",quote="",header=TRUE)
FuncID<-read.table("GOFuncIDs.cts_75.Mapping.padded.main.cumulative.summary_table.tsv",sep="\t",quote="",header=TRUE)
PFamID<-read.table("PFamIDs.cts_75.Mapping.padded.main.cumulative.summary_table.tsv",sep="\t",quote="",header=TRUE)
TIGRID<-read.table("TIGRFamIDs.cts_75.Mapping.padded.main.cumulative.summary_table.tsv",sep="\t",quote="",header=TRUE)
Phyllo<-read.table("Phyllo_Genus.cts_75.combined.main.summary_table.xls",sep="\t",quote="",header=T)
factors<-read.csv("phyllo_factors.csv",header=T)[-c(9:12),]
datasets<-list(ECID,ProcID,FuncID,PFamID,TIGRID,Phyllo)
remove(ECID,ProcID,FuncID,PFamID,TIGRID,Phyllo) #Dump old variables to save memory
#Grab max control values
CheckControls<-function(data){
  controls<-data[9:12,2:length(data)]
  A<-apply(controls,2,max)
  return(A)
}
#Clean data by subtracting max control value from each (also takes out control rows)
AdjustData<-function(data){
  data.max<-CheckControls(data)
  newdata<-data[-c(9:12),]
  for (i in 2:length(newdata)){
    newdata[i]<-newdata[i]-data.max[i-1]
    newdata[i][newdata[i]<0]<-0
    }
  return(newdata)
}
datasets<-lapply(datasets,AdjustData)


#Determine which observations correspond to which cities
CA.obs<-which(factors$city=="CA")
CA.group<-factors[CA.obs,]
HF.obs<-which(factors$city=="HF")
HF.group<-factors[HF.obs,]
DE.obs<-which(factors$city=="DE")
DE.group<-factors[DE.obs,]
Citylist<-list(CA.group,HF.group,DE.group)
for (i in 1:length(datasets)){
  datasets[[i]]$Total<-NULL
}
#Transpose and further manipulate data as needed
datasets<-lapply(datasets,t)
newdat<-list()
#Need to re-read the transposed datasets so that the sample IDs are recognized as headers
for (i in seq_along(datasets)){
  rownames(datasets[[i]])[1]<-""
  filename<-paste(i,".csv",sep="")
  write.table(datasets[[i]], file=paste("dataset",filename, sep="-"),sep=",",quote=TRUE,row.names=TRUE,col.names=FALSE)
  newdat[[i]]<-read.table(paste("dataset",filename,sep="-"),sep=",",header=T)
  unlink(paste("dataset",filename,sep="-"))
}
remove(datasets)

#Impute the average of the other drought replications to the third drought replications
for (i in 1:length(newdat)){
  newdat[[i]][,19]<-NA
  colnames(newdat[[i]])[19]<-"PHYLLO30"
  newdat[[i]]$PHYLLO30<-round((newdat[[i]]$PHYLLO28+newdat[[i]]$PHYLLO29)/2)
}

#Add new drought phenotype entry to HF group and factor file
HF.obs<-c(HF.obs,18)
HF.group[6,]<-factors[18,]<-NA
Samples<-c(as.character(factors$Sample_ID[1:17]),"PHYLLO30")
Samples.HF<-c(as.character(HF.group$Sample_ID[1:5]),"PHYLLO30")
factors$Sample_ID<-as.factor(Samples)
HF.group$Sample_ID<-as.factor(Samples.HF)
HF.group$treatment[6]<-factors$treatment[18]<-as.factor("drought")

#Subset data based on city (for some reason, indices are one off, but this is easily fixed)
Citysplitter<-function(data){
  list1<-list()
  list1[[length(list1)+1]]<-data[CA.obs+1]
  list1[[length(list1)+1]]<-data[DE.obs+1]
  list1[[length(list1)+1]]<-data[HF.obs+1]
  list1<-lapply(list1,function(y) cbind(data$X,y))
  return(list1)
}

uberlist<-lapply(newdat,Citysplitter)
#NOTE: Uberlist goes in the following order for cities: CA, DE, HF
row.names(CA.group)<-row.names(HF.group)<-row.names(DE.group)<-NULL
CA.group$city<-HF.group$city<-DE.group$city<-NULL
#Depending on the dataset you wish to run (either with subsets or without), you should run one of the following loops:
#This first loop is for the dataset with subsets based on city

for (i in 1:length(uberlist)){
  for (j in 1:length(uberlist[[1]])){
    filename<-paste(i,j,"csv",sep=".")
    #Set up group depending on which index j is:
    group<-switch(j, CA.group, DE.group, HF.group)
    #Some datasets have issues with normalization method, hence switch to method 2
    try(TwoStage_Package(uberlist[[i]][[j]],group,paste("sigtest","method2",filename,sep="-"),2))
    try(TwoStage_Package(uberlist[[i]][[j]],group,paste("sigtest",filename,sep="-"),1))
    print(c(i, j)) #Print the "coordinate" of the dataset so we know which ones specifically produc1e errors or insignificant results
  }
}

#This second loop is for analyzing the entire dataset, without accounting for the city
setwd("~/Metagenomic Data/Full Outputs") #Move outputs to new folder
rownames(factors)<-factors$city<-NULL
for (k in 1:length(newdat)){
  try(TwoStage_Package(newdat[[k]],factors,paste("sigtest",k,"method2","csv",sep="."),2))
  try(TwoStage_Package(newdat[[k]],factors,paste("sigtest",k,"csv",sep="."),1))
  print(k)
}