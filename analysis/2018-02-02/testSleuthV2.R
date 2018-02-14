## now we want to test diffex
source("../../bin/RNASeqDiffEx.R")

# --------------------------------------------------------------
## first get total data frame
# --------------------------------------------------------------
full.df <- buildDiffExDf()

# --------------------------------------------------------------
# test file corruption or extract valid HDF5 files 
# the link below explainsthe file corruption issue.
# https://support.hdfgroup.org/HDF5/release/known_problems/previssues.html
# --------------------------------------------------------------
test.file <- lapply(seq_along(full.df$path) , function(i){
  tryCatch(rhdf5::h5read(full.df$path[[i]],"aux/ids"), error = function(e) NULL)
})
file.position <- which(unlist(lapply(test.file, function(x) is.null(x))) %in% TRUE)
full.df[file.position,]


# --------------------------------------------------------------
# just want to compare hets to pos
# --------------------------------------------------------------
pos.hets <- subset(full.df, nf1Genotype != '-/-')
# error explained when we include sex inc = c(sex, 'Culture') as a covariate with this small sample set 
# https://github.com/pachterlab/sleuth/issues/69
fullmod <- buildSleuthModel(pos.hets, inc = c('Culture'), test = 'nf1Genotype', alt = '+/-')

# --------------------------------------------------------------
# plot
# --------------------------------------------------------------
# fground.tf <- unique(tab[which(as.numeric(tab[,'qval'])<0.01),1])

# test must match buildSleuthModel() test param 
# check fullmod$tests
plotVals(fullmod, qval = 0.01, test = 'nf1Genotype+/+', alt = 'Primary', collapseByGene = TRUE, prefix = 'culture')
plotVals(fullmod, qval = 0.01, ttype = c('protein_coding'), test = 'nf1Genotype+/+', alt = 'Primary', collapseByGene = TRUE, prefix = 'culture')
# ----------------------------------------------------------------------
## then do one-allele test again.
# ----------------------------------------------------------------------
oamod <- buildSleuthModel(full.df, inc = c('sex','Culture'), test = 'OneAllele', alt = '+')

#tab<-getSleuthTable(fullmod,test='OneAllele',alt='+')
##wirte list for go
#write(sapply(tab[,'transcript'],function(x) unlist(strsplit(x,split='.',fixed=T))[1]),file='pvalSortedensTrans.txt')
#pcoding<-tab[grep('protein_coding',tab[,1]),]

# ----------------------------------------------------------------------
# plot 
# ----------------------------------------------------------------------
# test must match buildSleuthModel() test param 
# check oamod$tests
plotVals(oamod, qval = 0.01, test = 'OneAllele+', alt = '+', collapseByGene = TRUE, prefix = 'oneAllele')
plotVals(oamod, qval = 0.01, ttype = c('protein_coding'), test = 'OneAllele+', alt = '+', collapseByGene = TRUE, prefix = 'oneAllele')


# --------------------------------------------------------------
# (not found) pmod ? 
# --------------------------------------------------------------
# ptab <- getSleuthTable(pmod, test = 'OneAllele+', alt = '+')
# plotVals(pmod, qval = 0.1, prefix = 'primaryOnly')
# plotVals(pmod, qval = 0.1, ttype = c('protein_coding'), prefix = 'primaryOnly')
# 
# plotVals(pmod, qval = 0.01, prefix = 'primaryOnly')
# plotVals(pmod, qval = 0.01, ttype = c('protein_coding'), prefix = 'primaryOnly')

# --------------------------------------------------------------
# brainstorm region
# --------------------------------------------------------------
# tab<-getSleuthTable(fullmod,'Culture','primary')

##wirte list for go
#write(sapply(tab[,'transcript'],function(x) unlist(strsplit(x,split='.',fixed=T))[1]),file='pvalSortedensTrans.txt')

#pcoding<-tab[grep('protein_coding',tab[,1]),]

#allgenes=unique(tab[,'gene'])##pretty sure this preserves order
#write(allgenes,file='prim_immor_allgenes_ranked.txt')

##now let's get protein coding
#pcgenes<-unique(tab[grep('protein_coding',tab[,'target_id']),'gene'])
#write(pcgenes,file='prim_immor_all_protcod_genes_ranked.txt')
# --------------------------------------------------------------
# brainstorm region
# --------------------------------------------------------------