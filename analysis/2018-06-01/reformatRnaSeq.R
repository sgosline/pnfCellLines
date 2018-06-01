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

full.tab<-all.dat%>%rename(id=`RNASeq Data Gencode`)%>%inner_join(all.tpm,by="id")
#now join into two tidied datasets

this.script='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2018-06-01/reformatRnaSeq.R'

gene.ex.tab<-full.tab%>%select(Sample=`Sample Name`,Gene,Value=tpm)%>%unique()
write.table(gene.ex.tab,file='updatedNTAPRnaSeqTpmTidied.tsv',sep='\t',row.names=F,col.names=T)

nf.pheno.tab<-full.tab%>%select(Sample=`Sample Name`,Value=`Sample Genotype`)%>%unique()
nf.pheno.tab$Phenotype=rep('NF1 Genotype',nrow(nf.pheno.tab))
write.table(nf.pheno.tab,file='NF1GenotypeNtapCellLines.tsv',sep='\t',row.names=F,col.names=T)

parentId='syn8282028'
synStore(File('updatedNTAPRnaSeqTpmTidied.tsv',parentId=parentId),used=all.dat$`RNASeq Data Gencode`,executed=this.script)

synStore(File('NF1GenotypeNtapCellLines.tsv',parentId=parentId),used=list('syn8397154'),executed=this.script)
