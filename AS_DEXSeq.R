#########################################
###   Load the required libraries     ###
#########################################

# DEXSeq
if(!require("DEXSeq", quietly = TRUE)){
  BiocManager::install("DEXSeq")
  library("DEXSeq")
}

# dplyr
if(!require("dplyr", quietly = TRUE)){
  install.packages("dplyr",dependencies=TRUE)
  library("dplyr")
}

# Multi cores
if(!require("BiocParallel", quietly = TRUE)){
  BiocManager::install("BiocParallel")
  library("BiocParallel")
}
BPPARAM = MulticoreParam(workers=2)

source("~/Desktop/Birmingham/Module6/Project_1/code/R/Subread_to_DEXSeq-master/load_SubreadOutput.R")

#####################################
###   Input counts into DEXSeq   ###
####################################

# Set the working directory
path <- "/Users/zhifan/Downloads/heart"
setwd(path)

# Input the featurecounts data matrix and sample information
samp <- data.frame(row.names = c("kidney_ERR315383", "kidney_ERR315443", 
                                 "kidney_ERR315468", "kidney_ERR315494",
                                 "heart1_ERR315328", "heart2_ERR315389", 
                                 "heart3_ERR315435"),
                   condition = c(rep("kidney",4), rep("heart",3))) # rep(c("kidney", "brain"), each=3)

dxd.fc <- DEXSeqDataSetFromFeatureCounts("kid_heart.count.txt", 
                                         flattenedfile = "hg38.ensGene.flat.gtf",
                                         sampleData = samp)

# Subset to 6 TSPAN genes of interest
genesForSubset <- read.table("geneSubset.txt", row.names = 1, stringsAsFactors=FALSE)[[1]]

dxd <- dxd.fc[geneIDs( dxd.fc ) %in% genesForSubset,]

# Check the data
colData(dxd)
head(counts(dxd),10)
head( featureCounts(dxd), 5 )

#############################################
###  Testing for differential exon usage  ###
############################################

# Normalisation
dxd <- estimateSizeFactors( dxd )

# Dispersion estimation
dxd <- estimateDispersions(dxd, fitType='local', BPPARAM=BPPARAM)

# plots the per-exon dispersion estimates versus the mean normalised count
plotDispEsts(dxd)

# Testing for differential exon usage
dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
res <- DEXSeqResults(dxd)

# Check the results
mcols(res)$description

# plots
plotDEXSeq( res, "ENSG00000134198", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( res, "ENSG00000134198", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( res, "ENSG00000134198", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( res, "ENSG00000134198", expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# Output a HTML file for visualization
DEXSeqHTML( res, genes = genesForSubset, color=c("#FF000080", "#0000FF80") )

# Clean the output table and save it
logname <- grep(names(res), pattern = 'log2fold', value = TRUE)
res.clean <- as(res[, c('groupID', 'featureID', 'exonBaseMean', 
                        logname, 'dispersion', 'stat', 'pvalue', 'padj')], 'data.frame')
names(res.clean)<- c("EnsemblID", "exonID", "meanBase", "log2FoldChange", 
                     "dispersion", "stat", "pvalue")
res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr') 
res.clean$chromosome <- as.character(seqnames( res$genomicData))
res.clean$exon.start <- start(res$genomicData)
res.clean$exon.end <- end(res$genomicData)
res.clean <- res.clean[, c("EnsemblID", "exonID", "meanBase", 
                           "log2FoldChange", "dispersion", "stat", "pvalue", 
                           "FDR", "chromosome", "exon.start", "exon.end")]  ### reorder the names nicely
res.clean <- res.clean[ order(res.clean$pvalue),]
# Store the output
write.csv(res.clean, quote = FALSE,
          file="DEXSeq_result.txt",
          row.names = FALSE)

sessionInfo()