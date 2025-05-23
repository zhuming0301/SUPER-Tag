setwd("projects/spatial_ct/")
set.seed(8)
suppressMessages(library(ArchR))
suppressMessages(library(parallel))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(EnsDb.Mmusculus.v79))
addArchRGenome("mm10")
addArchRThreads(threads = 24)
source("Code/ArchRtoSignac.R")
outputDirectory = "Save/MergedObject"

# 
ArrowFiles <- c(
    "E10.5" = "Save/E10.5/ArrowFiles/E10.5.arrow",
    "E11.5" = "Save/E11.5/ArrowFiles/E11.5.arrow",
    "E12.5" = "Save/E12.5/ArrowFiles/E12.5.arrow",
    "E13.5" = "Save/E13.5/ArrowFiles/E13.5.arrow",
    "E14.5" = "Save/E14.5/ArrowFiles/E14.5.arrow"
    )
ArchRProj <- ArchRProject(ArrowFiles = ArrowFiles, 
        outputDirectory = outputDirectory,
        copyArrows = F, showLogo = FALSE
        )
saveArchRProject(ArchRProj, outputDirectory = outputDirectory)
ArchRProj <- loadArchRProject(outputDirectory)

ArchRProj <- addIterativeLSI(
    ArchRProj, 
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    dimsToUse = 1:30,
    clusterParams = list(
        resolution = c(1.2), 
        maxClusters = 20, 
        sampleCells = 10000, 
        n.start = 10), 
    varFeatures = 25000, 
    force = TRUE
    )
ArchRProj <- addHarmony(
    ArchRProj = ArchRProj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)
ArchRProj <- addUMAP(ArchRProj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
    )
ArchRProj <- addUMAP(
    ArchRProj = ArchRProj, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", 
    force = TRUE
    )
ArchRProj <- addClusters(
    input = ArchRProj, 
    reducedDims = "IterativeLSI", 
    method = "Seurat", 
    name = "archr_clusters", 
    resolution = 1.2, 
    maxClusters = 20,
    force = TRUE
    )
ArchRProj <- addClusters(
    input = ArchRProj, 
    reducedDims = "Harmony", 
    method = "Seurat", 
    name = "harmony_clusters", 
    resolution = 1.2, 
    maxClusters = 20,
    force = TRUE
    )
p1 <- plotGroups(ArchRProj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin", 
    alpha = 0.4, 
    addBoxPlot = TRUE
    )
p2 <- plotGroups(ArchRProj,
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin", 
    alpha = 0.4, 
    addBoxPlot = TRUE
    )
p3 <- plotEmbedding(ArchRProj = ArchRProj, 
    colorBy = "cellColData", 
    name = "archr_clusters",
    size = 2, 
    embedding = "UMAP"
    )
p4 <- plotEmbedding(ArchRProj = ArchRProj, 
    colorBy = "cellColData", 
    name = "Sample", 
    size = 2, 
    embedding = "UMAP"
    )
p31 <- plotEmbedding(
    ArchRProj, 
    colorBy = "cellColData", 
    name = "Sample", 
    size = 2, 
    embedding = "UMAP_Harmony"
    )
p41 <- plotEmbedding(
    ArchRProj, 
    colorBy = "cellColData", 
    name = "harmony_clusters", 
    size = 2, 
    embedding = "UMAP_Harmony"
    )
p5 <- plotFragmentSizes(ArchRProj)
p5$labels$x <- "CUT&Tag Fragment Size (bp)"
p6 <- plotTSSEnrichment(ArchRProj)

plotPDF(p1, p2, p3, p4, p31, p41, p5, p6,
    name = "ArchR_QC_plots.pdf", 
    ArchRProj = ArchRProj, 
    addDOC = FALSE, 
    width = 4, 
    height = 4
    )
# %% Call peaks
ArchRProj <- addGroupCoverages(
    ArchRProj, 
    groupBy = "harmony_clusters", 
    maxFragments = 1*10^8, 
    force = TRUE
    )
ArchRProj <- addReproduciblePeakSet(
    ArchRProj,
    groupBy = "harmony_clusters",     
    pathToMacs2 = "macs2", 
    force = TRUE
    )
peaks <- getPeakSet(ArchRProj)
saveRDS(peaks, file = "../Save/ArchR_peaks.rds")
rtracklayer::export(
    peaks, 
    con = "Save/ArchR_peaks.bed", 
    format = "bed"
    )
# %% Cluster on PeakMatrix
# ArchRProj <- addPeakSet(ArchRProj, peaks, force = T)
ArchRProj <- addPeakMatrix(
    ArchRProj, 
    ceiling = 100, 
    binarize = FALSE,
    force = TRUE
    )
ArchRProj <- addIterativeLSI(
    ArchRProj = ArchRProj, 
    useMatrix = "PeakMatrix", 
    name = "IterativeLSI",  
    iterations = 2, 
    clusterParams = list(
        resolution = c(1.2), 
        sampleCells = 10000, 
        maxClusters = 20,
        n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:30, 
    force = TRUE
    )
ArchRProj <- addHarmony(
    ArchRProj = ArchRProj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)
ArchRProj <- addUMAP(
    ArchRProj = ArchRProj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", 
    force = TRUE
    )
ArchRProj <- addUMAP(
    ArchRProj = ArchRProj, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", 
    force = TRUE
    )
ArchRProj <- addClusters(
    input = ArchRProj, 
    reducedDims = "IterativeLSI", 
    method = "Seurat", 
    name = "archr_clusters", 
    resolution = 1.2, 
    maxClusters = 50,
    force = TRUE
    )
ArchRProj <- addClusters(
    input = ArchRProj, 
    reducedDims = "Harmony", 
    method = "Seurat", 
    name = "harmony_clusters", 
    resolution = 1.2, 
    maxClusters = 20,
    force = TRUE
    )
p7 <- plotEmbedding(
    ArchRProj, 
    colorBy = "cellColData", 
    name = "Sample", 
    size = 2, 
    embedding = "UMAP"
    )
p8 <- plotEmbedding(
    ArchRProj, 
    colorBy = "cellColData", 
    name = "archr_clusters", 
    size = 2, 
    embedding = "UMAP"
    )
p9 <- plotEmbedding(
    ArchRProj, 
    colorBy = "cellColData", 
    name = "Sample", 
    size = 2, 
    embedding = "UMAP_Harmony"
    )
p10 <- plotEmbedding(
    ArchRProj, 
    colorBy = "cellColData", 
    name = "harmony_clusters", 
    size = 2, 
    embedding = "UMAP_Harmony"
    )
plotPDF(p7, p8, p9, p10, 
    name = "ArchR_Peak_UMAP.pdf", 
    ArchRProj = ArchRProj, 
    addDOC = FALSE, 
    width = 6, 
    height = 6
    )
ArchRProj <- addImputeWeights(ArchRProj)
ArchRProj <- addCoAccessibility(
    ArchRProj, 
    reducedDims = "IterativeLSI"
    )

saveArchRProject(ArchRProj, outputDirectory = outputDirectory)

# To Seurat Object
ArchRProj <- loadArchRProject(outputDirectory)
SeuratObject <- ArchRtoSignacSpatial(
        ArchRProj, 
        fragments.files = c(
            "E10.5" = "Data/E10.5/fragments.tsv.gz",
            "E11.5" = "Data/E11.5/fragments.tsv.gz",
            "E12.5" = "Data/E12.5/fragments.tsv.gz",
            "E13.5" = "Data/E13.5/fragments.tsv.gz",
            "E14.5" = "Data/E14.5/fragments.tsv.gz"
            ), 
        spatial.files.dirs = c(
            "E10.5" = "Data/E10.5/spatial",
            "E11.5" = "Data/E11.5/spatial",
            "E12.5" = "Data/E12.5/spatial",
            "E13.5" = "Data/E13.5/spatial",
            "E14.5" = "Data/E14.5/spatial"
        ), 
        genome = "mm10"
        )
SeuratObject <- addDim(ArchRProject, SeuratObject, addUMAPs = "UMAP_Harmony", reducedDims = "Harmony")
SeuratObject <- addGeneScoreAssay(ArchRProject, SeuratObject)
cluster <- SeuratObject$harmony_clusters
SeuratObject$harmony_clusters <- factor(cluster, 
    levels = unique(cluster[order(as.numeric(sub("C", "", cluster)))]))
saveRDS(SeuratObject, file = "Save/MergedObject.rds")

DefaultAssay(SeuratObject) <- "genescore"
SeuratObject[['peaks']] <- NULL
SeuratObject@reductions <- SeuratObject@reductions["umap"]
SeuratObject[["genescore"]]$count <- NULL
saveRDS(SeuratObject, file = "Save/MergedObject.genescore.rds")

#  Add scRNA-seq annotation
MOCA_E11 <- readRDS("/home/disk/zhuming/projects/spatial_ct/public/sci/MOCA_E11.obj.rds")

ArchRProj <- addGeneIntegrationMatrix(
    ArchRProj = ArchRProj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = MOCA_E11,
    addToArrow = TRUE,
    groupRNA = "Main_cell_type",
    force = TRUE
)
getAvailableMatrices(ArchRProj)
getPeakSet(ArchRProj)
ArchRProj <- addPeak2GeneLinks(
    ArchRProj = ArchRProj,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneIntegrationMatrix"
)

saveArchRProject(ArchRProj)
