##build random forest model

library(randomForest)
library(synapser)
library(tidyverse)
require(nplr)
synLogin()

drug.dat<-subset(read.csv(synGet('syn17462699')$path,header=T),study_synapse_id=='syn4939906')
geneex.dat<-subset(read.csv(synGet('syn18421359')$path,header=T),study=='pNF Cell Line Characterization')

drug.response<-function(drug.dat){
  
  
}

findVariableDrugs<-function(drug.dat){
  res=drug.dat%>%group_by(DT_explorer_internal_id)%>%mutate(variance=var(response))%>%select(DT_explorer_internal_id,variance)%>%unique()%>%arrange(desc(variance))
  subset(res,variance!=Inf)
}

buildTestModel<-function(geneex.dat,drug.cats){
  require(randomForest)
  
}

binDrugResponse<-function(drug.dat){
  
  
}