##test survival analysis code in TCGA using gene expression signature

source("../../../RASPathwaySig/bin/sigSurvAnalysis.R",chdir=TRUE)
library(synapseClient)
library(tidyverse)
library(pheatmap)

this.script<-'https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2017-12-05/nf1SurvAnalysis.R'
synapseLogin()

geneVals<-read.csv(synGet("syn11274043")@filePath)
geneVals$NFStatus<-rep("Up in NF1 Mutants",nrow(geneVals))
geneVals$NFStatus[which(geneVals$b<0)]<-"Down in NF1 Mutants"

upgenes<-as.character(unique(select(subset(geneVals,b>0),gene))[,1])
downgenes<-as.character(unique(select(subset(geneVals,b<0),gene))[,1])


##now loop through up/down genes to create figures
for(dis in tcga.cancer.types){
  mut.sig<-survivalAnalysisByMutation(dis,upgenes,listname='Upregulated in NF1 mutants')
  expr.sig<-survivalAnalysisByExpression(dis,upgenes,listname='Upregulated in NF1 mutants')
  print(paste('Significance of upreg genes in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
  mut.sig<-survivalAnalysisByMutation(dis,downgenes,listname='Downreg in NF1 mutants')
  expr.sig<-survivalAnalysisByExpression(dis,downgenes,listname='Downreg in NF1 mutants')
  print(paste('Significance of downreg genes in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
  
  combined=c(upgenes,downgenes)
  
  mut.sig<-survivalAnalysisByMutation(dis,combined,listname='Diffex NF1 mutants')
  expr.sig<-survivalAnalysisByExpression(dis,combined,listname='Diffex in NF1 mutants')
  print(paste('Significance of NF1 genes in',dis,'is',expr.sig$pval,'(expression) and',mut.sig$pval,'(mutation)'))
}