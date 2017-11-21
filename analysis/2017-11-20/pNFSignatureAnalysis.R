##collect pNF signature and compare it to external datasets

library(synapseClient)
library(tidyverse)
library(pheatmap)

this.script<-''
synapseLogin()

geneVals<-read.csv(synGet("syn11274043")@filePath)
geneVals$NFStatus<-rep("Up in NF1 Mutants",nrow(geneVals))
geneVals$NFStatus[which(geneVals$b<0)]<-"Down in NF1 Mutants"

#get gene annotations
nfstat<-select(geneVals,gene,NFStatus)%>%unique
nfstat<-nfstat$NFStatus%>%setNames(nfstat$gene)

cnfData<-read.table(synGet("syn5714951")@filePath)

#get patient annotations
res<-synTableQuery('select distinct id,specimenID,individualID from syn9702734')@values
res<-res[which(res$id%in%colnames(cnfData)),]

patNames<-res$individualID%>%setNames(res$id)


##plot in cNF data
pheatmap(cnfData[intersect(rownames(cnfData),as.character(genes)),],annotation_row=data.frame(NFStatus=nfstat),annotation_col = data.frame(Patient=patNames))
ggsave('pNFGenesIncNFdata.png')

synStore(File('pNFGenesIncNFdata.png',parentId=''),executed=list(this.script),used=list(list('syn11274043'),list('syn5714951')))
##plot in TCGA data 


##plot in CCLE data


