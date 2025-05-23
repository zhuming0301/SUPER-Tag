# ArchRtoSignac
.getPeakMatrix <- function (ArchRProject) 
{
    print("In Progress:")
    print("Get Matrix From ArchRProject")
    peak_matrix <- ArchR::getMatrixFromProject(ArchRProject, 
        useMatrix = "PeakMatrix", binarize = FALSE)
    pm <- assays(peak_matrix)$PeakMatrix
    rownames(pm) <- paste0(as.character(seqnames(ArchRProject@peakSet)), 
        "-", as.character(start(ArchRProject@peakSet)), "-", 
        as.character(end(ArchRProject@peakSet)))
    print("Return peak matrix")
    pm
}
addDim <- function (ArchRProject, SeuratObject, addUMAPs = "UMAP", reducedDims = "IterativeLSI") {
    print("In Progress:")
    print("add UMAP From ArchRProject to SeuratObject")
    umap_df <- ArchRProject@embeddings[[addUMAPs]]$df %>% as.matrix
    dim(umap_df)
    colnames(umap_df) <- c("UMAP_1", "UMAP_2")
    SeuratObject@reductions$umap <- Seurat::CreateDimReducObject(embeddings = umap_df, 
        assay = "peaks")
    print("In Progress:")
    print("add reduction From ArchRProject to SeuratObject")
    if (reducedDims == "Harmony") {
        harmony_matrix <- ArchRProject@reducedDims$Harmony$matDR
        colnames(harmony_matrix) <- paste0("LSI_", 1:ncol(harmony_matrix))
        SeuratObject@reductions$harmony <- Seurat::CreateDimReducObject(embeddings = harmony_matrix, 
            assay = "peaks")
    }
    else if (reducedDims == "IterativeLSI") {
        LSI_matrix <- ArchRProject@reducedDims$IterativeLSI$matSVD
        colnames(LSI_matrix) <- paste0("LSI_", 1:ncol(LSI_matrix))
        SeuratObject@reductions$IterativeLSI <- Seurat::CreateDimReducObject(embeddings = LSI_matrix, 
            assay = "peaks")
    }
    else if (reducedDims == "IterativeLSI2") {
        LSI_matrix2 <- ArchRProject@reducedDims$IterativeLSI2$matSVD
        colnames(LSI_matrix2) <- paste0("LSI_", 1:ncol(LSI_matrix2))
        SeuratObject@reductions$IterativeLSI2 <- Seurat::CreateDimReducObject(embeddings = LSI_matrix2, 
            assay = "peaks")
    }
    print("Return SeuratObject")
    SeuratObject
}

addGeneScoreAssay <- function (ArchRProject, SeuratObject) {
    print("In Progress:")
    print("Get Gene Score Matrix From ArchRProject")
    GeneScore_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix = "GeneScoreMatrix")
    gsm <- assays(GeneScore_matrix)$GeneScoreMatrix
    print("Get Impute Weights Matrix From ArchRProject")
    ImputeWeights <- getImputeWeights(ArchRProject)
    impute_Matrix <- imputeMatrix(
      mat = gsm,
      imputeWeights = ImputeWeights,
      threads = 1,
      verbose = FALSE,
      logFile = createLogFile("imputeMatrix")
    )
    print("get Gene Features From ArchRProject")
    row.names(impute_Matrix) = getFeatures(ArchRProject)
    print("Return Gene Score Matrix")
    SeuratObject[['genescore']] <- CreateAssayObject(counts = impute_Matrix)
    SeuratObject
}
addGeneIntegrationAssay <- function (ArchRProject, SeuratObject) {
    print("In Progress:")
    print("Get Gene Integration Matrix From ArchRProject")
    GeneIntegration_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix = "GeneIntegrationMatrix")
    gim <- assays(GeneIntegration_matrix)$GeneIntegrationMatrix
    print("Get Impute Weights Matrix From ArchRProject")
    ImputeWeights <- getImputeWeights(ArchRProject)
    impute_Matrix <- imputeMatrix(
      mat = gim,
      imputeWeights = ImputeWeights,
      threads = 1,
      verbose = FALSE,
      logFile = createLogFile("imputeMatrix")
    )
    row.names(impute_Matrix) = elementMetadata(GeneIntegration_matrix)$name
    print("Return Gene Integration Matrix")
    SeuratObject[['RNA']] <- CreateAssayObject(counts = impute_Matrix)
    SeuratObject
}

.getAnnotations <- function(genome){
    if (genome == "mm10") {
        annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        seqlevelsStyle(annotations) <- "UCSC"
        genome(annotations) <- "mm10"
    }
    return(annotations)
}
ArchRtoSignacSpatial <- function(
  ArchRProject, 
  fragments.files, 
  spatial.files.dirs, 
  genome = "mm10", 
  annotation = NULL
  ) {
    samples <- unique(ArchRProject@cellColData$Sample)
    fragments.files.dir <- lapply(fragments.files, function(x) paste0(dirname(x), "/"))
    if(is.null(annotation)){annotation <- .getAnnotations(genome)}
    pm <- .getPeakMatrix(ArchRProject)
    print("Prepare Seurat list for each sample")
        seurat_list <- lapply(samples, function(sample) {
          print(paste0("Creating Seurat Object for ", sample))
          cur_pm <- pm[, grepl(paste0("^", sample, "#"), colnames(pm))]
          colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), "#"))[, 2]
          cur_meta <- ArchRProject@cellColData %>% as.data.frame %>% 
              subset(Sample == sample)
          rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), "#"))[, 2]
          cur_chromatin <- Signac::CreateChromatinAssay(counts = cur_pm, 
              sep = c("-", "-"), fragments = fragments.files[[sample]], 
              ranges = ArchRProject@peakSet, genome = genome, 
              annotation = annotation)
          cur_atac <- Seurat::CreateSeuratObject(cur_chromatin, 
              assay = "peaks", meta.data = cur_meta)
          image <- Seurat::Read10X_Image(
              image.dir = spatial.files.dirs[[sample]], 
              filter.matrix = TRUE
              )
          image <- subset(image, cells = sub(".*#", "", getCellNames(ArchRProject)))
          DefaultAssay(image) <- 'Spatial'
          cur_atac[[sample]] <- image
          cur_atac <- RenameCells(cur_atac, new.names = paste0(sample, "#", Cells(cur_atac)))
        })

    SeuratObject <- if (length(seurat_list) > 1) {
        print("Merging Objects...")
        merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
        }
        else {
            seurat_list[[1]]
        }
    print("Return SeuratObject")
    print(SeuratObject)
    SeuratObject
}

addMotifAssay <- function (ArchRProject, SeuratObject) {
    print("In Progress:")
    print("Get Motif Matrix From ArchRProject")
    motif_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix = "MotifMatrix")
    mtfm <- assays(motif_matrix)$z
    print("Get Impute Weights Matrix From ArchRProject")
    ImputeWeights <- getImputeWeights(ArchRProject)
    impute_Matrix <- imputeMatrix(
      mat = mtfm,
      imputeWeights = ImputeWeights,
      threads = 1,
      verbose = FALSE,
      logFile = createLogFile("imputeMatrix")
    )
    print("Return Motif Matrix")
    SeuratObject[['motif']] <- CreateAssayObject(counts = impute_Matrix)
    SeuratObject
}