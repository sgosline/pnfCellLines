##Today start toe valuate SNP segments using ABSOLUTE

#will it work with the pre-segmented data?
source("../../bin/CNVData.R")

segs<-tier1_segmentedData()
colnames(segs)<-c("ID",'Chromosome','Start','End','Num_Probes','Segment_Mean')

annots<-cnv_annotation_data()


library(ABSOLUTE)


run_abs<-function(fname){
    sigma.p <- 0
    max.sigma.h <- 0.02
    min.ploidy <- 0.90
    max.ploidy <- 8
    max.as.seg.count <- 5000
    max.non.clonal <- .4
    max.neg.genome <- .4
    copy_num_type <- "total"

    oname=gsub('.seg','_PNF_CellLineSample',fname)
    sname=paste(annots[match(gsub('.seg','',fname),annots$sample),c(5,3)],collapse=' ')

    RunAbsolute(fname,sample.name=sname,,min.ploidy=min.ploidy,max.ploidy=max.ploidy,
                max.sigma.h=max.sigma.h,platform='SNP_6.0',copy_num_type='total',
                min.mut.af=0,sigma.p=sigma.p,results.dir='abs_res',verbose=TRUE,
                output.fn.base=oname,primary.disease='PNF',max.as.seg.count=max.as.seg.count,
                max.non.clonal = max.non.clonal, max.neg.genome=max.neg.genome)

}

require(parallel)
mclapply(unique(segs$ID),function(x){
    fname=paste(x,'seg',sep='.')
    write.table(segs[which(segs$ID==x),-1],file=paste(x,'seg',sep='.'),sep='\t')
    run_abs(fname)
},mc.cores=6)
