##compare in vivo screens to pNF data
source("../../bin/ncatsSingleAgentScreens.R")
require(dplyr)
aucs<-getValueForAllCells('TAUC')
aucs<-as.data.frame(aucs)
aucs$Phenotype<-rownames(aucs)
auc.vals<-tidyr::gather(aucs,Sample,Response,1:(ncol(aucs)-1))%>%rename(TAUC=Response)%>%mutate(PhenSamp=paste(Phenotype,Sample,sep='_'))%>%select(TAUC,PhenSamp)

lac50<-getValueForAllCells('LAC50')
lac50<-as.data.frame(lac50)
lac50$Phenotype<-rownames(lac50)
lac.vals<-tidyr::gather(lac50,Sample,Response,1:(ncol(lac50)-1))%>%rename(LAC50=Response)%>%mutate(PhenSamp=paste(Phenotype,Sample,sep='_'))%>%select(LAC50,PhenSamp)

maxr<-getValueForAllCells("MAXR")
maxr<-as.data.frame(maxr)
maxr$Phenotype<-rownames(maxr)
maxr.vals<-tidyr::gather(maxr,Sample,Response,1:(ncol(maxr)-1))%>%rename(MAXR=Response)%>%mutate(PhenSamp=paste(Phenotype,Sample,sep='_'))%>%select(MAXR,PhenSamp)

##now get data
drug.effic<-synTableQuery('SELECT NCATSname,Efficacy FROM syn8339862')@values
drug.effic<-rename(drug.effic,Drug=NCATSname)

##now do all the joins into a single table so we can plot
full.tab<-inner_join(auc.vals,lac.vals,by='PhenSamp')%>%inner_join(maxr.vals,by='PhenSamp')%>%tidyr::separate(PhenSamp,into=c("Drug","Sample"),sep='_')

tab.with.effic<-left_join(full.tab,drug.effic,by='Drug')

##do one last tidy
fin.tab<-tidyr::gather(tab.with.effic,key='Measurement',value='Value',TAUC,LAC50,MAXR)

##NOW we can plot it! huzzah!!
library(ggplot2)

par.id<-'syn8339995'
for(m in unique(fin.tab$Measurement)){
  ggplot(filter(fin.tab,!is.na(Efficacy),Measurement==m)%>%filter(Sample%in%c("ipNF95.6","ipNF05.5 (single clone)","ipNF05.5 (mixed clone)")))+geom_boxplot(aes(x=Efficacy,y=Value,col=Sample))
  fname=paste('inVivoComparisonOf95.6And05.5CellLinesOnly',m,'.png',sep='')
  ggsave(fname)
  synStore(File(fname,parentId=par.id),executed=list(list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/bin/ncatsSingleAgentScreens.R'),list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2017-03-23/compareToInvivo.R')))

}
write.table(fin.tab,file='filteredCellLineResponseData.tsv',sep='\t',row.names=F,col.names=T,quote=F)
synStore(File('filteredCellLineResponseData.tsv',parentId=par.id),executed=list(list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/bin/ncatsSingleAgentScreens.R'),list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2017-03-23/compareToInvivo.R')))



#+facet_grid(Measurement~.)

