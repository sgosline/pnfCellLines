##run random forest
require(synapser)
synLogin()
synId='syn17462699'
require(tidyverse)
tab<-read.csv(synGet(synId)$path)%>%dplyr::rename(internal_id='DT_explorer_internal_id')

drug.map<-synTableQuery('SELECT distinct internal_id,std_name FROM syn17090819')$asDataFrame()

tab.with.id<-tab%>%left_join(drug.map,by='internal_id')

all.compounds<-unique(tab.with.id$std_name)
all.models<-unique(tab.with.id$model_name)
print(paste('Loaded',length(all.compounds),'compound response data over',length(all.models),'models'))

##now get gene data


source("../../bin/randomForestExpressionModel.R")
