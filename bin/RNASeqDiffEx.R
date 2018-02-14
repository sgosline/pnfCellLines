##now we want to test diffex
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("tidyverse")
usePackage("pheatmap")
usePackage("ggplot2")
#usePackage("reshape2")
#usePackage("dplyr")


# source("http://bioconductor.org/biocLite.R")
# biocLite("devtools")    # only if devtools not yet installed
# biocLite("pachterlab/sleuth")

# source("https://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

library('sleuth')
library('synapseClient')
library('rhdf5')
library(ggplot2)
library(pheatmap)
library(tidyverse)
require(reshape2)
require(dplyr)

synapseLogin()
# -------------------------------------------------------------
##this function will download and construct the data frame used to build
##the sleuth object
##includes all the variables needed
#' @param outliers - list of cell line names to expcloe
buildDiffExDf <- function(outliers=c()) {
  
  h5files <- synTableQuery("SELECT id,individualID,nf1Genotype,sex FROM syn7817226 where parentId='syn5579785' and fileFormat='h5'")@values
  # h5files <- allfiles[grep('*h5',allfiles[,1]),]
  
  ol <- which(h5files$individualID %in% outliers)
  if (length(ol) > 0)
    h5files = h5files[-ol,]
  
  ##build the contrast matrix from the annotations
  h5files$path <- sapply(h5files$id, function(id) synGet(id)@filePath)
  # h5files
  
  #df=data.frame(sample=unlist(snames),path=unlist(filepaths),Sex=factor(sex),Origin=factor(unlist(sors),levels=rev(unique(unlist(sors)))),Genotype=factor(unlist(sgens),levels=c("++","+-","--")))
  one.all <- rep('+', nrow(h5files))
  one.all[which(h5files$nf1Genotype == '-/-')] <- '-'
  h5files$OneAllele <- factor(one.all, levels = c('-','+'))
  
  culture <- synTableQuery("SELECT distinct 'Sample Name','Sample Origin', 'Media' FROM syn8397154 where \"RNASeq Data\" is not null")@values
  
  df <- h5files %>%
    inner_join(culture, by = c('individualID' = 'Sample Name')) %>%
    rename(Culture = 'Media', Origin = 'Sample Origin', sample = 'individualID')
  
  df <- df[-which(df$id %in% df$id[duplicated(df$id)] & is.na(df$Culture)),]
  
  return(df)
}

buildGencodeTargetMap <- function(df){
  
  sf <- sleuth_prep(df, ~ OneAllele + Culture + Origin)
  
  ## now get gene names
  t2g <- t(sapply(unique(sf$obs_raw$target_id),function(x) {
    res = unlist(strsplit(x, split = '|',fixed = T))
    list(transcript = res[1], gene = res[6])}))
  t2g <- data.frame(target_id = unique(sf$obs_raw$target_id),t2g)
  
  t2g <- apply(t2g, 2, unlist)
  write.table(t2g,file = '../../data/gencodeGeneTranscriptMap.csv', sep = ',', col.names = T, row.names = F)
}

##now we need to build the sleuth object
buildSleuthModel <- function(df,inc=c('sex','Culture'), test = 'OneAllele', alt = '+', prefix = ''){
  
  df <- df[,c(inc, test, 'sample', 'path')]
  
  # which covariate contains missing values 
  # remove them from analysis 
  missing.vals <- colnames(df)[colSums(is.na(df)) > 0]
  missing.vals <- missing.vals[which(missing.vals %in% inc)]
  if (length(missing.vals) > 0) {
    # print(paste('removing covairate columns with missing value from analysis: ',  missing.vals))
    # df <- df[,which(!names(df) %in% missing.vals)]
    # inc <- inc[which(!inc %in% missing.vals)]
    print(paste('removing rows with NA from covairate columns: ',  missing.vals))
    df <- df[complete.cases(df), ]
  }
  
  t2g2 <- read.table('../../data/gencodeGeneTranscriptMap.csv', sep = ',', header = T)
  
  ##first build sleuth object
  formstring <- paste('~', test)
  ex <- paste(inc, collapse = ' + ')
  if (ex != "") {
    formstring <- paste(formstring, ex , sep = ' + ')
  }
  
  # test file corruption or extract valid HDF5 files
  test.file <- lapply(seq_along(df$path) , function(i){
    tryCatch(rhdf5::h5read(df$path[[i]],"aux/ids"), error = function(e) NULL)
  })
  
  # remove corrupted file from list 
  file.position <- which(unlist(lapply(test.file, function(x) is.null(x))) %in% TRUE)
  if (length(file.position) > 0) {
    print(paste('removing corrupt file at dataframe position: ', file.position))
    df <- df[!unlist(lapply(test.file, function(x) is.null(x))),]
  }else{
    print('all files were identified as valid HDF5 files')
  }
  
  print(paste('Building sleuth model for', formstring))
  sf <- sleuth_prep(df, as.formula(formstring), target_mapping = t2g2)
  
  # now fit the full model 
  ffit <- sleuth_fit(sf, as.formula(formstring))
  
  # now perform wald-test question
  # https://en.wikipedia.org/wiki/Wald_test (not recomended on single simple full models)
  tt <- paste(test, alt, sep = '')
  if (!paste(test, alt, sep = '') %in% colnames(ffit$design_matrix)) {
    print(ffit$design_matrix)
    # find the test string in model.matrix columns 
    tt <- colnames(ffit$design_matrix)[grepl(test, colnames(ffit$design_matrix))]
    print(paste('Matrix model:', tt))
    # This function computes the Wald test on one specific 'beta' coefficient on every transcript.
    fp <- sleuth_wt(ffit, tt)
  }else{
    print(paste('Matrix model:', tt))
    fp <- sleuth_wt(ffit, tt)
  }
  
  return(fp)
}

getSleuthTable <- function(sleuthfit,test,alt){
  # write results to table
  res.tab <- sleuth_results(sleuthfit,paste(test, alt, sep = ''))
  rtab <- apply(res.tab, 2, unlist)
  rtab
}

plotGenesInSamples <- function(obj, transcripts, units = "tpm",
                               genes = NULL, annotes = NULL, fname = 'selectedGenesInData.pdf',
                               test = 'OneAllele', alt='+', collapseByGene = FALSE){
  beta <- test
  jm <- obj$obs_norm
  jm$gene <- obj$target_mapping$gene[match(jm$target_id, obj$target_mapping$target_id)]
  
  tabd_df <- jm[jm$target_id %in% transcripts, ]
  
  # select gene unit of measure to plot
  if (units == "tpm") {
    tabd_df <- dplyr::select(tabd_df, gene,target_id, sample,
                             tpm)
    if (collapseByGene) {
      tabd_df <- reshape2::dcast(tabd_df, gene ~ sample,
                                 value.var = "tpm", fun.aggregate = sum)
      annotes <- NULL
    }
    else
      tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample,
                                 value.var = "tpm")
  }else if (units == "est_counts") {
    tabd_df <- dplyr::select(tabd_df, gene,target_id, sample,
                             est_counts)
    if (collapseByGene) {
      tabd_df <- reshape2::dcast(tabd_df, gene ~ sample,
                                 value.var = "est_counts", fun.aggregate = sum)
      annotes <- NULL
    }else{
      tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample,
                                 value.var = "est_counts")
    }
  }else{
    stop("Didn't recognize the following unit: ", units)
  }
  
  dm <- design_matrix(obj)
  cul <- rep('Immortalized', nrow(dm))
  names(cul) <- rownames(dm)
  
  if ('CulturePrimary' %in% colnames(dm)) {
    cul[which(dm[,'CulturePrimary'] == 1)] <- 'Primary'
  }else{
    cul <- NULL
  }
  
  if ('OneAllele+' %in% colnames(dm)) {
    oa <- rep("-", nrow(dm))
    names(oa) <- rownames(dm)
    oa[which(dm[, beta] == 1)] <- '+'
    adf <- data.frame(oneNF1 = oa)
  }else{
    genotype <- rep("++", nrow(dm))
    names(genotype) <- rownames(dm)
    genotype[which(dm[,beta] == 1)] <- '+-'
    genotype[which(!dm[,beta] == 1)] <- '--'
    adf <- data.frame(genotype = genotype)
  }
  
  if (!is.null(cul)) {
    adf$culture <- cul 
  }
  
  if (!is.null(genes) && !collapseByGene) {
    rownames(tabd_df) <- genes[tabd_df[,1]]
  }else{
    rownames(tabd_df) <- tabd_df[,1] 
  }
  
  
  if (!is.null(annotes)) {
    names(annotes) <- rownames(tabd_df) 
  }
  
  # save heatmaps to output file path 
  if (is.null(annotes)) {
    pheatmap(log2(tabd_df[,-1] + 0.01), cellheight = 10, cellwidth = 10,
             annotation_col = adf, clustering_distance_rows = 'correlation',
             clustering_distance_cols = 'correlation', filename = fname)
    dev.off()
  }else{
    pheatmap(log2(tabd_df[,-1] + 0.01), cellheight = 10, cellwidth = 10,
             annotation_col = adf, clustering_distance_rows = 'correlation',
             clustering_distance_cols = 'correlation',
             annotation_row = data.frame(Type = annotes), filename = fname)
    dev.off()
  }
  
}


plotSingleGene <- function(obj, gene, prefix = '') {
  jm <- obj$obs_norm
  jm$gene <- obj$target_mapping$gene[match(jm$target_id, obj$target_mapping$target_id)]
  
  
  tabd_df <- jm[which(jm$gene == gene), ]
  
  tsum <- daply(tabd_df,'sample', function(td){sum(td$tpm)})
  newdf <- data.frame(TPM = tsum, sample = names(tsum))
  newdf$Genotype <- obj$sample_to_covariates$Genotype[match(names(tsum),obj$sample_to_covariates$sample)]
  newdf$Culture <- obj$sample_to_covariates$Culture[match(names(tsum),obj$sample_to_covariates$sample)]
  newdf$Origin <- obj$sample_to_covariates$Origin[match(names(tsum),obj$sample_to_covariates$sample)]
  
  pdf(paste(prefix,gene,'expression.pdf', sep = ''))
  p <-  ggplot(newdf) + geom_bar(aes(Origin,TPM, fill = Culture), stat = 'identity', position = 'dodge')
  print(p)
  dev.off()
}

plotVals <- function(data.obj, qval, ttype = c(), prefix = '', test = 'OneAllele', alt = '+', collapseByGene = FALSE){
  
  # test param in this function has to match buildSleuthModel() test param.
  res <- sleuth_results(data.obj, test = test)
  sel = which(res$qval < qval)
  print(paste("Found", length(sel),'diff ex transcripts at q=', qval))
  
  targs <- as.character(res$target_id)#[sel])
  trans.type <-  sapply(as.character(targs), function(x) {
    arr <- unlist(strsplit(x,split = '|', fixed = T))
    arr[length(arr)]
  })
  
  if (length(ttype) > 0) {
    sel <- intersect(sel, which(trans.type %in% ttype))
    print(length(sel))
  }
  
  targs <- targs[sel]
  gene <- res$gene[sel]
  trans <- res$transcript[sel]
  tnames <- paste(gene, trans, sep = '_')
  
  names(tnames) <- targs
  names(trans.type) <- tnames
  
  plotGenesInSamples(data.obj, targs,'tpm', tnames, trans.type[sel], 
                     fname = paste(prefix,'diffex',paste(ttype,collapse = '_'), 
                                   ifelse(collapseByGene,'Genes','Transcripts'), 'Q', qval, '.png', sep = ''), 
                     test, alt, collapseByGene = collapseByGene)
  
}

# General analysis
makeTablesAndPlots <- function(model, test, alt, prefix = ''){
  fname <- paste(test,'variable', alt, 'test', prefix, sep = '_')
  ##1
  tab <- getSleuthTable(model, test, alt)
  #store table
  write.table(tab, file = paste(fname,'SleuthResults.csv',sep = ''), sep = ',', quote = F)
  sigs <- which(as.numeric(tab[,'qval']) < 0.1)
  
  pc <- sigs[grep('protein_coding',tab[sigs,'target_id'])]
  
  pcg <- unique(tab[pc,'gene'])
  print(paste("Found",length(sigs),'differentially expressed transcripts, of which',length(pc),'are protein coding representing',length(pcg),'unique genes'))
  
  pdf(paste(fname,'VolcanoPlot.pdf', sep = ''))
  res <- plot_volcano(model, paste(test, alt, sep = ''), point_alpha = 0.8)
  print(res)
  dev.off()
  
  pdf(paste(fname,'PCAPlot.pdf', sep = ''))
  res <- plot_pca(model, text_labels = TRUE, color_by = test, point_size = 8)
  print(res)
  dev.off()
  
  gvals <- tab[,'gene']
  names(gvals) <- tab[,'target_id']
  plotGenesInSamples(model, tab[pc,'target_id'], units = "tpm",
                     genes = NULL, annotes = NULL, collapseByGene = TRUE,
                     fname = paste(fname,'sigGenesByTpmInHeatmap.pdf'), test = test, alt = alt)
}

getGOList <- function(sleuthfit){
  
}