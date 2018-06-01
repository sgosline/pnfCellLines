##reformat tidied expression counts of all pNF cell culture cell lines with genotype

require(synapser)
synLogin()

require(tidyverse)
##first get all cell lines, genotype, culture conditions

all.dat<-synapser::synTableQuery('SELECT `Sample Name`,`Sample Genotype`,`Media`,`RNASeq Data Gencode` FROM syn8397154')$asDataFrame()%>%subset(!is.na(`RNASeq Data Gencode`))

##then pull TPM files

all.tpm<-do.call('rbind',lapply(all.dat$`RNASeq Data Gencode`,function(x){
  tab<-NULL
  try(tab<-read.table(synGet(x)$path,header=T)%>%select(Gene=HugoSymbol,tpm))
  if(is.null(tab))
    tab<-read.table(synGet(x)$path,header=T)%>%select(Gene=HGNCSymbol,tpm)
  new.tab<-tab%>%select(Gene,tpm)%>%group_by(Gene)%>%mutate(totalTpm=sum(tpm))%>%select(Gene,tpm=totalTpm)%>%unique()
  new.tab$id=rep(x,nrow(new.tab))
  return(new.tab)
  }
  ))

#now join into two tidied datasets
gene.ex.file<-''

nf.pheno.file<-''

parentId='syn8282028'

