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


fragments.files <- c(
    "E10.5" = "Data/E10.5/fragments.tsv.gz",
    "E11.5" = "Data/E11.5/fragments.tsv.gz",
    "E12.5" = "Data/E12.5/fragments.tsv.gz",
    "E13.5" = "Data/E13.5/fragments.tsv.gz",
    "E14.5" = "Data/E14.5/fragments.tsv.gz"
    )
spatial.files.dirs <- c(
    "E10.5" = "Data/E10.5/spatial",
    "E11.5" = "Data/E11.5/spatial",
    "E12.5" = "Data/E12.5/spatial",
    "E13.5" = "Data/E13.5/spatial",
    "E14.5" = "Data/E14.5/spatial"
)
fragments.files <- sapply(
    fragments.files,
    function(x) normalizePath(x, mustWork = TRUE))
spatial.files.dirs <- sapply(
    spatial.files.dirs,
    function(x) normalizePath(x, mustWork = TRUE))
# %%
if (!dir.exists('./ArchR/')) {dir.create('./ArchR/')}
if (!dir.exists('./Save/')) {dir.create('./Save/')}
for(sample in names(fragments.files)){
    image <- Seurat::Read10X_Image(
        image.dir = spatial.files.dirs[[sample]], 
        filter.matrix = TRUE
        )
    setwd("./ArchR")
    outputDirectory = paste0("../Save/", sample)

    ArrowFiles <- createArrowFiles(
        inputFiles = fragments.files[sample], 
        sampleNames = sample, 
        validBarcodes = Seurat::Cells(image),
        minTSS = 0, minFrags = 0, 
        maxFrags = 1e+07, force = FALSE
        )
    ArchRProj <- ArchRProject(ArrowFiles = ArrowFiles, 
        outputDirectory = outputDirectory,
        copyArrows = TRUE, showLogo = FALSE
        )
    # ArchRProj <- loadArchRProject(outputDirectory)
    ArchRProj <- addIterativeLSI(ArchRProj, 
        clusterParams = list(resolution = c(1.2), maxClusters = 20, sampleCells = 10000, n.start = 10),
        force = TRUE
        )
    ArchRProj <- addUMAP(ArchRProj, nNeighbors = 30, force = TRUE)
    ArchRProj <- addClusters(input = ArchRProj, resolution = 1.2, force = TRUE)
    p1 <- plotGroups(ArchRProj, groupBy = "Sample", colorBy = "cellColData", 
        name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
    p2 <- plotGroups(ArchRProj, groupBy = "Sample", colorBy = "cellColData",
        name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
    p3 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "Clusters",size = 2, embedding = "UMAP")
    p4 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "Sample", size = 2, embedding = "UMAP")

    ArchRProj <- addGroupCoverages(ArchRProj, maxFragments = 1*10^8, force = TRUE)
    ArchRProj <- addReproduciblePeakSet(ArchRProj, pathToMacs2 = "macs2", force = TRUE)
    peaks <- getPeakSet(ArchRProj)
    saveRDS(peaks, file = paste0("../Save/ArchR_peaks.", sample, ".rds"))
    rtracklayer::export(peaks, con = paste0("../Save/ArchR_peaks.", sample, ".bed"), format = "bed")
    ArchRProj <- addPeakMatrix(ArchRProj, ceiling = 100)
    ArchRProj <- addIterativeLSI(ArchRProj = ArchRProj, useMatrix = "PeakMatrix", force = TRUE,
        clusterParams = list(resolution = c(1.2), maxClusters = 20, sampleCells = 10000, n.start = 10))
    ArchRProj <- addUMAP(ArchRProj, nNeighbors = 30, force = TRUE)
    ArchRProj <- addClusters(input = ArchRProj, resolution = 1.2, force = TRUE)
    
    p5 <- plotFragmentSizes(ArchRProj)
    p5$labels$x <- "CUT&Tag Fragment Size (bp)"
    p6 <- plotTSSEnrichment(ArchRProj)
    p7 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "Sample", size = 2, embedding = "UMAP")
    p8 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "Clusters", size = 2, embedding = "UMAP")
    plotPDF(p1, p2, p3, p4, p5, p6, p7, p8, 
        name = "ArchR_Plots.pdf", 
        ArchRProj = ArchRProj, addDOC = FALSE, width = 4, height = 4)
    ArchRProj <- addImputeWeights(ArchRProj)
    ArchRProj <- addCoAccessibility(ArchRProj, reducedDims = "IterativeLSI")
    setwd("../")
    saveArchRProject(ArchRProj, paste0("Save/", sample))
    metadata = getCellColData(ArchRProj)
    write.csv(metadata, file = paste0("Save/ArchRProj.cellcoldata.", sample, ".csv"), quote = F)
    getGroupBW(ArchRProj)
    }


MOCA_E11 <- readRDS("/home/disk/zhuming/projects/spatial_ct/public/sci/MOCA_E11.obj.rds")
for(sample in names(fragments.files)) {
    outputDirectory <- paste0("Save/", sample)
    ArchRProj <- loadArchRProject(outputDirectory)

    SeuratObject <- ArchRtoSignacSpatial(
        ArchRProj, 
        fragments.files, 
        spatial.files.dirs, 
        genome = "mm10", 
        annotations
        )
    SeuratObject <- addDim(ArchRProject, SeuratObject, addUMAPs = "UMAP", reducedDims = "IterativeLSI")
    SeuratObject <- addGeneScoreAssay(ArchRProject, SeuratObject)
    # transfer anchors
    transfer.anchors <- FindTransferAnchors(reference = MOCA_E11, 
        query = SeuratObject, query.assay = 'genescore', reduction = 'cca', dims = 1:30)
    predicted.labels <- TransferData(anchorset = transfer.anchors,
    refdata = MOCA_E11$Main_cell_type, weight.reduction = SeuratObject[['IterativeLSI']], dims = 2:30)
    SeuratObject <- AddMetaData(object = SeuratObject, metadata = predicted.labels)
    cluster <- SeuratObject$Clusters
    SeuratObject$Clusters <- factor(cluster, 
        levels = unique(cluster[order(as.numeric(sub("C", "", cluster)))]))
    saveRDS(SeuratObject, file = paste0("Save/SeuratObject.", sample, ".rds"))
}
