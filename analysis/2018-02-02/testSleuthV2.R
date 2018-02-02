##now we want to test diffex
source("../../bin/RNASeqDiffEx.R")


##first get total data frame
full.df<-buildDiffExDf()

#just want to compare hets to pos
fullmod<-buildSleuthModel(subset(full.df,nf1Genotype!='-/-'),inc=c('Sex','Culture'),test='nf1Genotype',alt='+/-')

#tab<-getSleuthTable(fullmod,'Culture','primary')

##wirte list for go
#write(sapply(tab[,'transcript'],function(x) unlist(strsplit(x,split='.',fixed=T))[1]),file='pvalSortedensTrans.txt')

#pcoding<-tab[grep('protein_coding',tab[,1]),]



#allgenes=unique(tab[,'gene'])##pretty sure this preserves order
#write(allgenes,file='prim_immor_allgenes_ranked.txt')

##now let's get protein coding
#pcgenes<-unique(tab[grep('protein_coding',tab[,'target_id']),'gene'])
#write(pcgenes,file='prim_immor_all_protcod_genes_ranked.txt')


##now try to plot
#fground.tf=unique(tab[which(as.numeric(tab[,'qval'])<0.01),1])

plotVals(fullmod,qval=0.01,test='Culture',alt='primary',collapseByGene=TRUE,prefix='culture')
plotVals(fullmod,qval=0.01,ttype=c('protein_coding'),test = 'Culture',alt='primary',collapseByGene=TRUE,prefix='culture')

####
##then do one-allele test again.

oamod<-buildSleuthModel(full.df,test='OneAllele',alt='+')

##
#tab<-getSleuthTable(fullmod,test='OneAllele',alt='+')

##wirte list for go
#write(sapply(tab[,'transcript'],function(x) unlist(strsplit(x,split='.',fixed=T))[1]),file='pvalSortedensTrans.txt')


#pcoding<-tab[grep('protein_coding',tab[,1]),]

##now try to plot

plotVals(oamod,qval=0.01,test='OneAllele',alt='+',collapseByGene=TRUE,prefix='oneallele')
plotVals(oamod,qval=0.01,ttype=c('protein_coding'),test='OneAllele',alt='+',collapseByGene=TRUE,prefix='oneallele')


##
ptab<-getSleuthTable(pmod,test='OneAllele',alt='+')
plotVals(pmod,qval=0.1,prefix='primaryOnly')
plotVals(pmod,qval=0.1,ttype=c('protein_coding'),prefix='primaryOnly')

plotVals(pmod,qval=0.01,prefix='primaryOnly')
plotVals(pmod,qval=0.01,ttype=c('protein_coding'),prefix='primaryOnly')


