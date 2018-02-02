##compare kinase target data to drug response in pNF cell lines


library(xlsx)
require(tidyverse)

##now get the pNF cell line response data
source("../../bin/ncatsSingleAgentScreens.R")
source('../../bin/singleDrugAnalysis.R')

##store kinase output once,then comment out
#syn<-synStore(File('nbt.2017-S2.xls',parentId='syn11600912'))
this.script='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2017-12-08/testKinaseModel.R'
tab<-read.xlsx(file=synGet('syn11600914')@filePath,sheetName='Sheet1',startRow=3,check.names=F)[,1:179]

aucMat<-getRecalculatedAUCMatrix()

rownames(aucMat)<-sapply(rownames(aucMat),function(x) gsub(' ','.',x))

drugs<-intersect(rownames(aucMat),colnames(tab))

#write the cas numbers to send to robert
#write.table(colnames(tab)[-1],file='casNumbers.txt',quote=F,row.names=F,col.names=F)

fullfile<-read.table(synGet('syn5522643')@filePath,header=T,as.is=T,sep=',',quote='"',comment.char = '')
ncats2name<-fullfile[,c(2,37)]

#now store robert's lookup, then comment out
synStore(File('ncats_to_cas_map.txt',parentId='syn11600912'))
ncats2cas<-read.table(synGet('syn11600920')@filePath,header=T,quote='"')
merged<-ncats2name%>%rename(ncats="NCGC.SID")%>%inner_join(ncats2cas)

sub.merged<-subset(merged,cas%in%colnames(tab))

##now tidy the table with the ncats names
tidied.tab<-tab%>%rename(Target="compound CAS#:")%>%gather(key="cas",value="kinaseResponse",2:ncol(tab))%>%inner_join(sub.merged,by='cas')%>%select(Target,name,kinaseResponse,cas)


##get mean changes b/w genotypes
genotypes<-synTableQuery("SELECT 'Sample Name','Sample Genotype' FROM syn8397154")@values
genevalues<-genotypes$`Sample Genotype`%>%setNames(genotypes$`Sample Name`)

genevalues[which(genevalues=='+/-')]<-'+/+'
genevalues$ipn02.8<-'+/+'
res<-aucDifferentialResponse(aucMat,genevalues[colnames(aucMat)])
res$Drug<-rownames(res)

#now we have the full data frame! let's do elastic net

##negative numbers indiciate GREATER expression in NF1 -/- samples
##logfc is +/+ / -/-
full.dat<-tidied.tab%>%inner_join(select(res,logFC,name="Drug"),by='name')
library(glmnet)
nas<-which(is.na(full.dat$kinaseResponse))
if(length(nas)>0)
  full.dat<-full.dat[-nas,]

#try this

mean.resp<-acast(full.dat,name~Target,value.var="kinaseResponse",fun.aggregate=mean,fill=100)
resps<-select(full.dat,name,logFC)%>%unique
resps<-resps$logFC[match(rownames(mean.resp),resps$name)]
res<-glmnet(x=mean.resp,y=resps,alpha=0.5)

##now try to visualize results
lambs<-res$lambda
coef.mat<-lapply(lambs,function(x) {
  coeffs=coef(res,s=x)
  res<-data.frame(Value=as.numeric(coeffs))
  res$lambda=rep(x,nrow(res))
  res$Target=rownames(coeffs)
  return(res)})

coef.df<-data.frame(do.call(rbind,coef.mat))
require(ggplot2)
nz.coefs<-coef.df[which(coef.df$Value!=0),]

ggplot(nz.coefs)+geom_line(aes(x=lambda,y=Value,col=Target))
ggsave('allNonZeroCoefficients.png')
synStore(File('allNonZeroCoefficients.png',parentId='syn11600912'),used=list(this.script))

no.int<-subset(nz.coefs,Target!='(Intercept)')

ggplot(no.int)+geom_line(aes(x=lambda,y=Value,col=Target))
ggsave('allNonZeroCoefficientsNoIntercept.png')
synStore(File('allNonZeroCoefficientsNoIntercept.png',parentId='syn11600912'),used=list(this.script))
thresh<-subset(no.int,lambda>.02)

ggplot(thresh)+geom_line(aes(x=lambda,y=Value,col=Target))
ggsave('allNonZeroCoefficientsNoInterceptLambdaover0.02.png')
synStore(File('allNonZeroCoefficientsNoInterceptLambdaover0.02.png',parentId='syn11600912'),used=list(this.script))

thresh<-subset(no.int,lambda>.04)

ggplot(thresh)+geom_line(aes(x=lambda,y=Value,col=Target))
ggsave('allNonZeroCoefficientsNoInterceptLambdaover0.04.png')
synStore(File('allNonZeroCoefficientsNoInterceptLambdaover0.04.png',parentId='syn11600912'),used=list(this.script))

#save parameters
#how well does model predict?
#plot AUC

pdf('lambdavsCor.pdf')
plot(lambs,sapply(lambs,function(x) cor(resps,predict(res,newx=mean.resp,s=x))),xlab='Lambda',ylab='Correlation')
dev.off()

synStore(File('lambdavsCor.pdf',parentId='syn11600912'),used=list(this.script))
