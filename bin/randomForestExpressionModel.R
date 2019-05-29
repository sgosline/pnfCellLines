##build random forest model

library(randomForest)
library(synapser)
library(tidyverse)
require(nplr)
synLogin()

#get drug data
drug.dat<-subset(read.csv(synGet('syn17462699')$path,header=T),study_synapse_id=='syn4939906')%>%subset(response_type=='AUC_Simpson')%>%subset(organism_name=='human')


std <- synTableQuery("SELECT  DT_explorer_internal_id, name FROM syn18506944")$asDataFrame() %>%
  distinct()%>%select(-c(ROW_ID,ROW_VERSION))

drug.dat<-drug.dat%>%left_join(std,by='DT_explorer_internal_id')

##get gene exs
geneex.dat<-subset(read.csv(synGet('syn18421359')$path,header=T),study=='pNF Cell Line Characterization')

ex.mat<-reshape2::acast(geneex.dat,Symbol~specimenID,value.var="zScore",fun.aggregate=function(x) mean(x,na.rm=T))

#ex.mat<-geneex.dat%>%select(specimenID,zScore,Symbol)%>%group_by(specimenID,zScore,Symbol)%>%mutate(mean=mean(zScore))%>%select(-zScore)%>%spread(key=Symbol,value=specimenID)

drug.response<-function(drug.dat){


}

#let's focus on teh drugs that are most variable!
findVariableDrugs<-function(drug.dat){
  res=drug.dat%>%group_by(name,response_type)%>%mutate(variance=var(response))%>%select(name,response_type,variance)%>%unique()%>%arrange(desc(variance))
  res<-subset(res,variance!=Inf)
  res
}

buildTestModel<-function(geneex.dat,drug.cats){
  require(randomForest)

}

##bin drug response by upper and lower third
binDrugResponse<-function(drug.dat,drug){
  vals=subset(drug.dat,name==drug)%>%distinct()#%>%subset(response_type==res$response_type[[2]])%>%distinct()
  quants=quantile(vals$response,probs=c(0.333,0.666))
  bins=sapply(vals$response,function(x) ifelse(x<quants[1],'LOW',ifelse(x>quants[2],'HIGH',NA)))
  data.frame(vals$model_name,bin=bins)
}
