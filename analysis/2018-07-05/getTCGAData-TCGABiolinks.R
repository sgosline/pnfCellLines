library(synapser)
synLogin()
library(plyr)
library(tidyverse)
library(data.table)
library(TCGAbiolinks)
library(Biobase)
# 
# GDCprepare <- function(query,
# save = FALSE,
# save.filename,
# directory = "GDCdata",
# summarizedExperiment = TRUE,
# remove.files.prepared = FALSE,
# add.gistic2.mut = NULL,
# mut.pipeline = "mutect2",
# mutant_variant_classification = c("Frame_Shift_Del",
# "Frame_Shift_Ins",
# "Missense_Mutation",
# "Nonsense_Mutation",
# "Splice_Site",
# "In_Frame_Del",
# "In_Frame_Ins",
# "Translation_Start_Site",
# "Nonstop_Mutation")){
# 
# isServeOK()
# if(missing(query)) stop("Please set query parameter")
# if(any(duplicated(query$results[[1]]$cases)) & query$data.type != "Clinical data" & query$data.type !=  "Protein expression quantification") {
# dup <- query$results[[1]]$cases[duplicated(query$results[[1]]$cases)]
# dup <- query$results[[1]][query$results[[1]]$cases %in% dup,c("cases","experimental_strategy")]
# dup <- dup[order(dup$cases),]
# print(knitr::kable(dup))
# stop("There are samples duplicated. We will not be able to preapre it")
# }
# if(!save & remove.files.prepared) {
# stop("To remove the files, please set save to TRUE. Otherwise, the data will be lost")
# }
# # We save the files in project/source/data.category/data.type/file_id/file_name
# source <- ifelse(query$legacy,"legacy","harmonized")
# files <- file.path(query$results[[1]]$project, source,
# gsub(" ","_",query$results[[1]]$data_category),
# gsub(" ","_",query$results[[1]]$data_type),
# gsub(" ","_",query$results[[1]]$file_id),
# gsub(" ","_",query$results[[1]]$file_name))
# 
# files <- file.path(directory, files)
# 
# if(!all(file.exists(files))) stop(paste0("I couldn't find all the files from the query. ",
# "Please check if the directory parameter right or GDCdownload downloaded the samples."))
# 
# if(grepl("Transcriptome Profiling", query$data.category, ignore.case = TRUE)){
# data <- readTranscriptomeProfiling(files = files,
# data.type = ifelse(!is.na(query$data.type),  as.character(query$data.type),  unique(query$results[[1]]$data_type)),
# workflow.type = unique(query$results[[1]]$analysis_workflow_type),
# cases = query$results[[1]]$cases,
# summarizedExperiment)
# } else if(grepl("Copy Number Variation",query$data.category,ignore.case = TRUE)) {
# data <- readCopyNumberVariation(files, query$results[[1]]$cases)
# }  else if(grepl("DNA methylation",query$data.category, ignore.case = TRUE)) {
# data <- readDNAmethylation(files, query$results[[1]]$cases, summarizedExperiment, unique(query$platform))
# }  else if(grepl("Protein expression",query$data.category,ignore.case = TRUE)) {
# data <- readProteinExpression(files, query$results[[1]]$cases)
# }  else if(grepl("Simple Nucleotide Variation",query$data.category,ignore.case = TRUE)) {
# if(grepl("Masked Somatic Mutation",query$data.type,ignore.case = TRUE) | source == "legacy")
# suppressWarnings(data <- readSimpleNucleotideVariationMaf(files))
# }  else if(grepl("Clinical|Biospecimen", query$data.category, ignore.case = TRUE)){
# data <- readClinical(files, query$data.type, query$results[[1]]$cases)
# summarizedExperiment <- FALSE
# } else if (grepl("Gene expression",query$data.category,ignore.case = TRUE)) {
# if(query$data.type == "Gene expression quantification")
# data <- readGeneExpressionQuantification(files = files,
# cases = query$results[[1]]$cases,
# summarizedExperiment = summarizedExperiment,
# experimental.strategy = unique(query$results[[1]]$experimental_strategy))
# 
# if(query$data.type == "miRNA gene quantification")
# data <- readGeneExpressionQuantification(files = files,
# cases = query$results[[1]]$cases,
# summarizedExperiment = FALSE,
# experimental.strategy = unique(query$results[[1]]$experimental_strategy))
# if(query$data.type == "miRNA isoform quantification")
# data <- readmiRNAIsoformQuantification(files = files,
# cases = query$results[[1]]$cases)
# 
# if(query$data.type == "Isoform expression quantification")
# data <- readIsoformExpressionQuantification(files = files, cases = query$results[[1]]$cases)
# 
# }
# # Add data release to object
# if(summarizedExperiment & !is.data.frame(data)){
# metadata(data) <- list("data_release" = getGDCInfo()$data_release)
# }
# 
# 
# if((!is.null(add.gistic2.mut)) & summarizedExperiment) {
# message("=> Adding GISTIC2 and mutation information....")
# genes <- tolower(levels(EAGenes$Gene))
# if(!all(tolower(add.gistic2.mut) %in% genes)) message(paste("These genes were not found:\n",
# paste(add.gistic2.mut[! tolower(add.gistic2.mut) %in% genes],collapse = "\n=> ")))
# add.gistic2.mut <- add.gistic2.mut[tolower(add.gistic2.mut) %in% tolower(genes)]
# if(length(add.gistic2.mut) > 0){
# info <- colData(data)
# for(i in unlist(query$project)){
# info <- get.mut.gistc.information(info,
# i,
# add.gistic2.mut,
# mut.pipeline = mut.pipeline,
# mutant_variant_classification = mutant_variant_classification)
# }
# colData(data) <- info
# }
# }
# if("samples" %in% colnames(data)){
# if(any(duplicated(data$sample))) {
# message("Replicates found.")
# if(any(data$is_ffpe)) message("FFPE should be removed. You can do data with the following command:\ndata <- data[,!data$is_ffpe]")
# print(as.data.frame(colData(data)[data$sample %in% data$sample[duplicated(data$sample)],c("is_ffpe"),drop=F]))
# }
# }
# 
# 
# if(save){
# if(missing(save.filename) & !missing(query)) save.filename <- paste0(query$project,gsub(" ","_", query$data.category),gsub(" ","_",date()),".RData")
# message(paste0("Saving file:",save.filename))
# save(data, file = save.filename)
# message("File saved")
# 
# # save is true, due to the check in the beggining of the code
# if(remove.files.prepared){
# # removes files and empty directories
# remove.files.recursively(files)
# }
# }
# 
# return(data)
# }

pancan <- TCGAbiolinks:::getGDCprojects()$project_id
pancan <- pancan[grep("TCGA", pancan)]

##get all files that need to be retrieved 
query <- lapply(pancan, function(x){
  GDCquery(project = x,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM-UQ")
})

##get all data
lapply(query, function(x){
  GDCdownload(x, method = "api", files.per.chunk = 10)
})

data <- lapply(query, function(x){
  GDCprepare(x, summarizedExperiment = FALSE)
})
names(data) <- pancan
saveRDS(data, "TCGA_expression_data.rds")

nf1_mut_samples <- lapply(pancan, function(x){
  y <- gsub("TCGA-", "", x)
  maf <- GDCquery_Maf(y, pipelines = "muse") %>% 
    mutate(Tumor_Sample_Barcode = gsub("normalized_count_","",Tumor_Sample_Barcode) %>% 
             gsub("\\D-[[:alnum:]]+-[[:alnum:]]+-[[:alnum:]]+$", "",.))
  
  nf1_mut <- maf %>% 
    filter(Entrez_Gene_Id == "4763", !grepl("silent", 
                                            Variant_Classification, ignore.case = TRUE))
  nf1_mut_samples <- unique(nf1_mut$Tumor_Sample_Barcode) 
  
  nf1_wt <- maf %>% 
    filter(!Tumor_Sample_Barcode %in% nf1_mut_samples) 
  
  nf1_wt_samples <- unique(nf1_wt$Tumor_Sample_Barcode)
  
  list("mutant" = nf1_mut_samples, "wt" = nf1_wt_samples)
})

names(nf1_mut_samples) <- pancan

saveRDS(nf1_mut_samples, "TCGA_NF1_nonsilent_mutation_status.rds")
