##combine annotations for all genes

library(synapseClient)
synapseLogin()

id.map<-synTableQuery('SELECT "Sample Name","Exome-Seq Identifier","Sample Genotype" FROM syn5014742')@values

anno.files<-synQuery("select name,id from entity where parentId=='syn6086887'")
anno.files<-anno.files[grep('txt',anno.files$entity.name),]

anno.files$seqID=sapply(anno.files$entity.name,function(x) unlist(strsplit(x,split='_'))[1])

parse.annovar<-function(file){
  tab<-read.table(file,sep='\t',header=T,comment.char = '#')
#  library(ggplot2)
  
  ##remove any variants for which all 3 gene models agree on a functional prediction
  #any.mis=which(apply(tab[,c('KnownGeneExonFunction','EnsgeneExonFunction','RefgeneExonFunction')],1,function(x) any(x=="")||any(x=='synonymous SNV')))
  #func.tab<-tab[-any.mis,]
  func.tab<-tab
  print(paste("Returning",nrow(func.tab),'variants in',length(unique(func.tab$KnownGeneGeneName)),'genes'))
  func.tab
}

all.func.tabs<-lapply(anno.files$entity.id,function(x) parse.annovar(synGet(x)@filePath))

full.tab<-do.call('rbind',lapply(anno.files$entity.id))
names(all.func.tabs)<-id.map$`Sample Name`[match(anno.files$seqID,id.map$`Exome-Seq Identifier`)]

