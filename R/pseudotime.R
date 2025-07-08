setwd("projects/spatial_ct/")
# %% RNA
set.seed(8)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

project <- readRDS("Save/ST.E14.5.Forebrain.rds")
project.sub <- subset(project, Clusters %in% c("R6","R16","R2","R13"))
cds <- as.cell_data_set(project.sub)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = Cells(subset(project, Clusters %in% c("R6"))))

pseudotime <- pseudotime(cds, reduction_method = "UMAP")
project.sub$pseudotime <- pseudotime

p <- SpatialPlot(project.sub, features = "pseudotime", pt.size.factor = 2, stroke = 0, combine = F)
p1 <- p[[1]] + scale_fill_gradientn(colours = ArchR::ArchRPalettes$horizonExtra)

plot_genes <- c("Pax6", "Tbr1", "Eomes", "Neurod2")
plot_cds <- cds[rownames(cds) %in% plot_genes, colData(cds)$Clusters %in% c("R6","R16","R2","R13")]
p2 <- monocle3::plot_genes_in_pseudotime(plot_cds,  panel_order = c("Pax6", "Eomes", "Tbr1", "Neurod2"), label_by_short_name = F, ncol = 2)
ggsave(p1 + p2, filename = "Plots/pseudotime.svg", width = 8, height = 4)

library(ggplot2)
library(scales)
library(viridis)
library(dplyr)

plot_gene_pseodutime <- function(project, gene){
    df <- FetchData(project, vars = c("pseudotime", gene)) %>%
        rename(expr = gene)
    df <- df %>%
    mutate(pt_scaled = scales::rescale(pseudotime, to = c(0, 1)))
    ggplot(df, aes(x = pseudotime, y = expr, color = pseudotime)) +
        geom_point(size = 0.2, stroke = 0) +
        scale_color_gradientn(colors = c("darkblue", "blue", "purple", "pink", "orange", "yellow")) +
        geom_smooth(color = "black", se = FALSE, method = "loess", size = 0.2) +
        labs(x = "PseudoTime", y = gene) +
        theme_bw(base_size = 5) +
        theme(legend.position = "bottom",
            panel.background = element_blank(),
            plot.background = element_blank(), 
            panel.grid = element_blank(),
            axis.text = element_text(size = 3)) +  
        guides(color = guide_colorbar(title = NULL, barwidth = 1, barheight = 0.2))
}
p1 <- plot_gene_pseodutime(project.sub, "Pax6") + theme(legend.position = "none", axis.title.x = element_blank())
p2 <- plot_gene_pseodutime(project.sub, "Eomes") + theme(legend.position = "none", axis.title.x = element_blank())
p3 <- plot_gene_pseodutime(project.sub, "Tbr1") + theme(legend.position = "none", axis.title.x = element_blank())
p4 <- plot_gene_pseodutime(project.sub, "Neurod2")
ggsave((p1+p2)/(p3+p4), filename = "Plots/pseudotime2.svg", width = 2, height = 1.5)

# %% CUT&Tag
library(ArchR)
outputDirectory = "Save/E14.5_Forebrain_H3K27ac"
ArchRProj = loadArchRProject(outputDirectory)

monocle_Cortex <- getMonocleTrajectories(
    ArchRProj = ArchRProj,
    name = "Cortex",
    useGroups = c("C6","C11","C13"),
    groupBy = "Clusters",
    principalGroup = "C6",
    embedding = "UMAP",
    clusterParams = list(k = 50),
    seed = 1
)
ArchRProj <- addMonocleTrajectory(
    ArchRProj = ArchRProj,
    name = "Cortex",
    useGroups = c("C6","C11","C13"),
    groupBy = "Clusters",
    monocleCDS = monocle_Cortex,
    force = TRUE
)

metadata <- getCellColData(ArchRProj)
write.csv(metadata, file = "Save/ArchRProj.metadata.Trajectories.csv")

project.CT <- readRDS("Save/SeuratObject.E14.5_Forebrain_H3K27ac.rds")

project.CT@meta.data$pseudotime <- metadata[Cells(project.CT), ]$Cortex
project.ct.sub <- subset(project.CT, cells = row.names(metadata[!is.na(metadata$Cortex), ]))
p <- SpatialPlot(project.ct.sub, features = "pseudotime", pt.size.factor = 1.2, stroke = 0, combine = F)
p3 <- p[[1]] + scale_fill_gradientn(colours = ArchR::ArchRPalettes$horizonExtra)
ggsave(p3, filename = "Plots/pseudotime.CT.spatial.svg", width = 5, height = 5)

pl1 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "GeneScoreMatrix", 
    name = "Pax6", continuousSet = "horizonExtra")
pl2 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "GeneScoreMatrix", 
    name = "Eomes", continuousSet = "horizonExtra")
pl3 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "GeneScoreMatrix", 
    name = "Tbr1", continuousSet = "horizonExtra")
pl4 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "GeneScoreMatrix", 
    name = "Neurod2", continuousSet = "horizonExtra")

plotPDF(pl1[[2]], pl2[[2]], pl3[[2]], pl4[[2]], name = "pseudotime.CT", width = 5, height = 5, addDOC = F)

pl5 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "MotifMatrix", 
    name = "z:Pax6_614", continuousSet = "horizonExtra")
pl6 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "MotifMatrix", 
    name = "z:Eomes_768", continuousSet = "horizonExtra")
pl7 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "MotifMatrix", 
    name = "z:Tbr1_769", continuousSet = "horizonExtra")
pl8 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "MotifMatrix", 
    name = "z:Neurod2_69", continuousSet = "horizonExtra")
plotPDF(pl5[[2]], pl6[[2]], pl7[[2]], pl8[[2]], name = "pseudotime.Motif.CT", width = 5, height = 5, addDOC = F)




# Trajectory heatmap 
trajGSM <- getTrajectory(ArchRProj, name = "Cortex", useMatrix = "GeneScoreMatrix", log2Norm = TRUE, trajectoryLabel = "Clusters")
p5 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), labelTop = 80,
    colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGSM)$label)))
plotPDF(p5, name = "Plot-Cortex-Traj-Heatmaps.pdf", addDOC = FALSE, width = 6, height = 12)

# pseudotime Enhancers
CREs <- readRDS("Save/CREs.rds")

ArchRProj <- addFeatureMatrix(
  input = ArchRProj,
  features = CREs,
  matrixName = "CREMatrix",
  ceiling = 10^9,
  binarize = FALSE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addFeatureMatrix")
)

crematrix <- getMatrixFromProject(ArchRProj, useMatrix = "CREMatrix")
pm <- assays(crematrix)$CREMatrix
rownames(pm) <- names(CREs)
ImputeWeights <- getImputeWeights(ArchRProj)
impute_Matrix <- imputeMatrix(
      mat = pm,
      imputeWeights = ImputeWeights,
      threads = 1,
      verbose = FALSE,
      logFile = createLogFile("imputeMatrix")
    )
cre_df <- as.data.frame(t(impute_Matrix))
cre_df <- cre_df[getCellNames(ArchRProj),]
for(i in colnames(cre_df)){
    ArchRProj <- addCellColData(ArchRProj, data = cre_df[[i]], name = i, cells = getCellNames(ArchRProj), force = TRUE)
}
getCellColData(ArchRProj)
plotList = c()
for(i in names(CREs)){
    pl1 <- plotTrajectory(ArchRProj, trajectory = "Cortex", colorBy = "cellColData", 
        name = i, continuousSet = "horizonExtra")
    plotList[[i]] <- pl1[[2]]
}
plotPDF(plotList, name = "pseudotime.enhancers.CT", width = 5, height = 5, addDOC = F)
