#########################################
###   Load the required libraries     ###
#########################################

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  library("BiocManager")
}

if(!require("ggplot2", quietly = TRUE)){
  install.packages("ggplot2",dependencies=TRUE)
  library("ggplot2")
}

# Library for GTF file
if(!require("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)){
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
}

# Library for AS prediction
if(!require("SGSeq", quietly = TRUE)){
  BiocManager::install("SGSeq")
  library("SGSeq")
}

#########################################
###   Input Annotation and BAM file   ###
#########################################

# Set the working directory.
path <- "/Users/zhifan/Downloads/heart"
setwd(path)

# Setup the sample info
sample_name <- c("braina_ERR315432","brainb_ERR315477","brainc_ERR315455") 
file_bam <- c("braina_ERR315432.bam", "brainb_ERR315477.bam",
              "brainc_ERR315455.bam")

# Get the BAM file information
bam_df <- data.frame(sample_name, file_bam, stringsAsFactors = FALSE)
sample_info <- getBamInfo(bam_df)

# Transcript annotation file
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- keepSeqlevels(txdb, c('chr1','chr7','chr11','chr12'))
# transcript features are extracted from the TxDb object
txf <- convertToTxFeatures(txdb)

# Genomic coordinates of 6 Tetraspanins
TSPAN2 <- GRanges("chr1", IRanges(start=115048011, end=115089503),
                  strand="-")
CD9 <- GRanges("chr12", IRanges(start=6199715, end=6238271),
               strand="+") # TSPAN29
TSPAN12 <- GRanges('chr7', IRanges(start=120787320,end=120858402),
                   strand='-')
CD81 <- GRanges("chr11", IRanges(start=2376177, end=2397419),
                strand="+") # TSPAN28
CD151 <- GRanges('chr11', IRanges(start=832887 ,end=839831),
                 strand='+') # TSPAN24
CD63 <- GRanges('chr12', IRanges(start=55725323 ,end=55729707),
                strand='-') # TSPAN30
tspan <- GRangesList(TSPAN2=TSPAN2,CD9=CD9,
                     TSPAN12=TSPAN12,CD81=CD81,
                     CD151=CD151,CD63=CD63)
# Only features overlapping the genomic locus of the tspan of interest
txf_tspan <- txf[txf %over% tspan]

#############################################################
###  Splice graph analysis based on annotated transcripts ###
#############################################################

# Converts transcript features to splice graph features 
sgfc_tspan <- analyzeFeatures(sample_info, features = txf_tspan, cores = 2)
sgfc_tspan
colData(sgfc_tspan)
rowRanges(sgfc_tspan)
head(counts(sgfc_tspan))
head(FPKM(sgfc_tspan))

# Splice graph analysis based on annotated transcripts
# FPKMs for splice graph features can be visualized with function plotFeatures.
# TSPAN2,CD9,TSPAN12,CD81,CD151,CD63
for (i in 1:length(tspan)){
  name <- paste0('brain_plot/',names(tspan)[i],'_FPKM_splice.png')
  png(name, width = 800, height = 800)
  plotFeatures(sgfc_tspan, which = tspan[[i]])
  dev.off()
}

############################################################
###   Splice graph analysis based on de novo prediction  ###
###########################################################

for (i in 1:length(tspan)){
  
  sgfc_pred <- analyzeFeatures(sample_info, which = tspan[[i]], cores = 2)
  
  # Annotate the predicted feature
  sgfc_pred <- annotate(sgfc_pred, txf)
  
  # Plot the splice graph and FPKMs
  name <- paste0('brain_plot/',names(tspan)[i],'_FPKM_annotation.png')
  png(name, width = 800, height = 800)
  plotFeatures(sgfc_pred, geneID = 1, color_novel = "red")
  dev.off()
} 

# For the coverage plot, have to do it manually
sgfc_pred <- analyzeFeatures(sample_info, which = TSPAN2)
sgfc_pred <- analyzeFeatures(sample_info, which = CD9)
sgfc_pred <- analyzeFeatures(sample_info, which = TSPAN12)
sgfc_pred <- analyzeFeatures(sample_info, which = CD81)
sgfc_pred <- analyzeFeatures(sample_info, which = CD151)
sgfc_pred <- analyzeFeatures(sample_info, which = CD63)

# Annotate the predicted feature
sgfc_pred <- annotate(sgfc_pred, txf)
#head(rowRanges(sgfc_pred))

# Splice variant identification
sgvc_pred <- analyzeVariants(sgfc_pred)
mcols(sgvc_pred)

# Splice variant quantification
variantFreq(sgvc_pred)
plotVariants(sgvc_pred, eventID = 1, color_novel = "red")

# per-base read coverages and splice junction counts can be visualized with function plotCoverage.
par(mfrow = c(5, 1), mar = c(1, 3, 1, 1))
plotSpliceGraph(rowRanges(sgfc_pred), geneID = 1, toscale = "none", color_novel = "red")
for (j in 1:dim(sgfc_pred)[2]) {
  plotCoverage(sgfc_pred[, j], geneID = 1, toscale = "none")
}

sessionInfo()
