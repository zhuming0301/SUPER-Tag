####### This code is for using FigR for gene regulation analysis
####### https://github.com/liranmao/Spatial_multi_omics
library(Seurat)
library(Signac)
library(ArchR)
library(SummarizedExperiment)
library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ComplexHeatmap)
library(networkD3)

setwd("projects/spatial_ct/")
# Preparation
# Extract counts matrix and features from ATAC
ArchRProj <- loadArchRProject("Save/MergedObject")
Peak_matrix <- getMatrixFromProject(ArchRProj, useMatrix='PeakMatrix')
pkm <- assays(Peak_matrix)$PeakMatrix
# Create SummarizedExperiment object
ATAC.se <- SummarizedExperiment(assays=list(counts=pkm), rowRanges = getPeakSet(ArchRProj))
dim(ATAC.se) # Peaks x Cells

# # RNA
Gene_matrix <- getMatrixFromProject(ArchRProj, useMatrix='GeneIntegrationMatrix')
gsm <- assays(Gene_matrix)$GeneIntegrationMatrix
row.names(gsm) = Gene_matrix@elementMetadata@listData$name
dim(gsm)


# ## FigR
# # Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                                 RNAmat = gsm,
                                 genome = "mm10", # One of hg19, mm10 or hg38 
                                 nCores = 10,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = 100)
saveRDS(cisCorr,file = 'Save/cisCorr_FigR.GeneIntegrationMatrix.rds')

cisCorr <- readRDS('Save/cisCorr_FigR.GeneIntegrationMatrix.rds')
cisCorr = cisCorr[grep("*Rik|^Gm|^mt-|^Rps|^Rpl|^Mir", cisCorr$Gene, invert = T),]


cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= 0.05)
pdf(paste0("Plots/DOCR.GeneIntegrationMatrix.pdf"), height = 4, width = 4)
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                    cutoff = 5, # No. sig peaks needed to be called a DORC, cutoff=2
                    labelTop = 25,
                    returnGeneList = TRUE, # Set this to FALSE for just the plot
                    labelSize = 2,
                    force = 10
                    )
dev.off()

markerPeaks <- readRDS("Save/markerpeaks.rds")

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

df.list = c()
genes = c()
for (i in seq_along(markerList)){
    list <- markerList[[i]]
    if (nrow(list)>0){
    peaks <- paste(list$seqnames, paste(list$start, list$end, sep = '-'), sep = ':')
    cisCorr.filt <- cisCorr[cisCorr$PeakRanges %in% peaks,] %>% filter(pvalZ <= 0.05)
    df <- as.data.frame(table(cisCorr.filt$Gene))
    df_sorted <- df[order(-df$Freq), ]
    df_sorted <- df_sorted[df_sorted$Freq >= 5,]
    df_sorted$cluster = paste0("C", i)
    df_sorted$gene = as.character(df_sorted$Var1)
    df.list[[i]] = df_sorted
    genes = c(genes, df_sorted$gene)
    print(head(df.list[[i]][df.list[[i]]$gene %in% uniqe_gene,], n=20))
    pdf(paste0("Plots/DOCR.cluster.C", i, ".pdf"), height = 4, width = 4)
    dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                        cutoff = 5, # No. sig peaks needed to be called a DORC, cutoff=2
                        labelTop = 25,
                        returnGeneList = TRUE, # Set this to FALSE for just the plot
                        labelSize = 2,
                        force=2)
    dev.off()
    }
}


