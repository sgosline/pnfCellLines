##run random forest
require(synapser)
require(tidyverse)

synLogin()

#get drug data


std <- synTableQuery("SELECT  DT_explorer_internal_id, name FROM syn18506944")$asDataFrame() %>%
  distinct()%>%select(-c(ROW_ID,ROW_VERSION))
drug.dat<-subset(read.csv(synGet('syn17462699')$path,header=T),study_synapse_id=='syn4939906')%>%subset(response_type=='AUC_Simpson')%>%subset(organism_name=='human')

drug.dat<-drug.dat%>%left_join(std,by='DT_explorer_internal_id')
drug.dat$model_name=gsub('b C','bC',drug.dat$model_name)
drug.dat$model_name=gsub('ipnNF95.11C','ipnNF95.11c',drug.dat$model_name)
drug.dat$model_name=gsub("ipNF05.5$","ipNF05.5 (single clone)",drug.dat$model_name)

##get gene exs
geneex.dat<-subset(read.csv(synGet('syn18421359')$path,header=T),study=='pNF Cell Line Characterization')

ex.mat<-reshape2::acast(geneex.dat,specimenID~Symbol,value.var="zScore",fun.aggregate=function(x) mean(x,na.rm=T))


source("../../bin/randomForestExpressionModel.R")

var.drugs<-findVariableDrugs(drug.dat)

require(pROC)
require(parallel)
rocs<-do.call(rbind,mclapply(unique(var.drugs$name),function(x){
    cats<-binDrugResponse(drug.dat,drug=x)
    auc.val=NA
    prval=NA
    if(is.null(cats))
        return(NA)
    mod<-NULL
    try(mod<-buildTestModel(ex.mat,cats))
    if(!is.null(mod)){
        try(res<-evalModel(mod,ex.mat,cats))
        auc.val=res$perf
        prval=res$pred

    }

    return(list(drug=x,numSamps=length(intersect(cats$vals.model_name,rownames(ex.mat))),AUC=auc.val,pred=prval))))
},mc.cores=30))

#now plot results
