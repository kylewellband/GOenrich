#!/usr/bin/env Rscript

args <- commandArgs(T)

if (length(args) != 3) {
    message(paste0("USAGE: RunTopGO.R <geneid2GO mapping file> <reference gene set> <test gene set>"))
    stop("You are missing an input file.")
}

if( suppressMessages(!require(topGO)) ) {
    if (!requireNamespace("BiocManager"))
        install.packages("BiocManager")
    BiocManager::install("topGO")
    suppressMessages(library(topGO))
}

gene2GO_file <- args[1]
refset_file <- args[2]
testset_file <- args[3]

testset <- read.table(testset_file, stringsAsFactors = F)[,1]
refset <- read.table(refset_file, stringsAsFactors = F)[,1]
gene2GO <- readMappings(gene2GO_file)

testset_w_go <- testset[testset %in% names(gene2GO)]
if ( length(testset_w_go) == 0 ) {
    stop("No genes in the test set have GO annotation! Double check your input files.")
}

refset_w_go <- refset[refset %in% names(gene2GO)]
if ( length(refset_w_go) == 0 ) {
    stop("No genes in the ref set have GO annotation! Double check your input files.")
}


gene2GO <- gene2GO[names(gene2GO) %in% refset_w_go]
gene_list <- factor(as.integer(refset_w_go %in% testset_w_go))
names(gene_list) <- refset_w_go

GOdata <- new("topGOdata",
              description = testset_w_go,
              ontology = "BP",
              allGenes = gene_list,
              nodeSize = 5,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

result_Fisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

write.table(GenTable(GOdata, result_Fisher, topNodes = length(result_Fisher@score)), 
            paste0("02results/", sub("(.*)?\\..*$", "\\1", basename(args[3])), "_GOenrich.txt"), 
            quote = F, row.names = F, sep = "\t")

#printGraph(GOdata, result_Fisher, firstSigNodes = 5, fn.prefix = paste0(base_dir, args[4]), pdfSW = T, useInfo = "np")

