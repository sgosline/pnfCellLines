##build random forest model

library(randomForest)
library(synapser)
library(tidyverse)
require(nplr)

#ex.mat<-geneex.dat%>%select(specimenID,zScore,Symbol)%>%group_by(specimenID,zScore,Symbol)%>%mutate(mean=mean(zScore))%>%select(-zScore)%>%spread(key=Symbol,value=specimenID)

drug.response<-function(drug.dat){


}

#let's focus on teh drugs that are most variable!
findVariableDrugs<-function(drug.dat){
  res=drug.dat%>%group_by(name,response_type)%>%mutate(variance=var(response))%>%select(name,response_type,variance)%>%unique()%>%arrange(desc(variance))
  res<-subset(res,variance!=Inf)
  res
}


evalModel<-function(mod,mat,cats){
    require(ROCR)
    res=NA
    pred=NULL
    try(pred<-predict(mod,mat,type='prob'))
    if(!is.null(pred)){
        samps=intersect(cats$vals.model_name,rownames(pred))
        print(samps)
        high.vals=cats[match(samps,cats$vals.model_name),2]=='HIGH'
        high.preds=pred[samps,'HIGH']
        pre=prediction(high.preds,high.vals)
        per=performance(pre,measure='auc')
        res=list(pred=pre,perf=per@y.values)
     }
    return(res)
}

buildTestModel<-function(geneex.dat,drug.cats){
    require(randomForest)
    drug.cats<-subset(drug.cats,!is.na(bin))
    samps<-intersect(drug.cats$vals.model_name,rownames(geneex.dat))
    if(length(samps)<2)
        return(NULL)
    print(samps)
    randomForest(x=geneex.dat[samps,],y=drug.cats$bin[match(samps,drug.cats$vals.model_name)])
}

##bin drug response by upper and lower third
binDrugResponse<-function(drug.dat,drug){
    vals=subset(drug.dat,name==drug)%>%distinct()%>%group_by(model_name)%>%mutate(meanresponse=mean(response,na.rm=T))%>%select(model_name,meanresponse)%>%distinct()
                                        #%>%subset(response_type==res$response_type[[2]])%>%distinct()
  quants=quantile(vals$meanresponse,probs=c(0.333,0.666))
  bins=sapply(vals$meanresponse,function(x) ifelse(x<quants[1],'LOW',ifelse(x>quants[2],'HIGH',NA)))
  data.frame(vals$model_name,bin=bins)
}
