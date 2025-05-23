setwd("projects/spatial_ct/")
source("Code/SpatialPlot_new.R")
options(scipen = 999)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
cols = ArchR::ArchRPalettes$stallion
names(cols) = paste0("C", 1:length(cols))
suppressMessages(extrafont::loadfonts())
if (!dir.exists('./Plots/')) {dir.create('./Plots/')}

DimPlot_theme <- theme(
    text = element_text(family = "Arial", size = unit(3, "pt")), 
    plot.title = element_text(hjust = 0.5, family = "Arial", size = unit(4, "pt"), face = "bold"),
    axis.title = element_text(family = "Arial", size = unit(3, "pt")), 
    axis.text = element_text(family = "Arial", size = unit(3, "pt")),  
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    legend.text = element_text(family = "Arial", size = unit(3, "pt")),
    legend.key.size = unit(4, 'pt'),
    aspect.ratio = 1
  )

Spatial_theme <- theme_void() + theme(
    text = element_text(family = "Arial", size = unit(4, "pt")),
    plot.title = element_text(hjust = 0.5, family = "Arial", size = unit(5, "pt"), face = "bold"),
    legend.text = element_text(family = "Arial", size = unit(4, "pt")), 
    legend.key.size = unit(4, 'pt'),
    aspect.ratio = 1
  )
# Single Projects
SeuratObjects = list()
metadata = data.frame()
for (i in c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5")) {
    SeuratObjects[[i]] <- readRDS(paste0("Save/SeuratObject.", i, ".rds"))
    df = SeuratObjects[[i]]@meta.data
    metadata = rbind(metadata, df)
}
write.csv(metadata, file = "Save/SeuratMetaData.csv", quote = F)

for (sample in names(spatial.files.dirs)){
    image <- Seurat::Read10X_Image(
        image.dir = spatial.files.dirs[[sample]], 
        filter.matrix = TRUE
        )
    image <- RenameCells(image, new.names = paste0(sample, "#", Cells(image)))
    DefaultAssay(image) <- 'Spatial'
    SeuratObjects[[sample]][[sample]]  <- image
}

plots1 = list()
plots2 = list()
plots3 = list()

for (i in 1:5) {
plots1[[i]] <- SpatialPlot(SeuratObjects[[i]], group.by = "Clusters", cols = cols,
  pt.size.factor = 1, stroke = 0, crop = F) + NoLegend()

plots2[[i]] <- DimPlot(SeuratObjects[[i]], group.by = "Clusters", cols = cols,
  pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes() + DimPlot_theme + ggtitle("")

plots3[[i]] <- DimPlot(SeuratObjects[[i]], group.by = "predicted.id", 
  pt.size = 0.1) + theme(aspect.ratio = 1, legend.position = 'none') + NoAxes() + DimPlot_theme + ggtitle("")

}
p1 <- cowplot::plot_grid(plotlist = plots1, ncol = 1)
p2 <- cowplot::plot_grid(plotlist = plots2, ncol = 1)
p3 <- cowplot::plot_grid(plotlist = plots3, ncol = 1)
combined_plot <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
ggsave(combined_plot, file = "Plots/UMAP-Clusters-2.svg", width = 6, height = 12)

# Merge Projects
metadata = read.csv("Save/SeuratMetaData.csv", row.names = 1)
MergeObject <- readRDS("Save/MergedObject.genescore.rds")
head(metadata)
# QC
plot <- VlnPlot(MergeObject,
      features = c("nFrags", "nCount_peaks", "FRIP", "TSSEnrichment", "BlacklistRatio"),
      group.by = "Sample", pt.size =0, combine = F)
    
plot[[1]] <- plot[[1]] + ylim(0,300000)
plot[[2]] <- plot[[2]] + ylim(0,150000)
plot <- lapply(plot, function(p) {p + 
    geom_boxplot(width = 0.2, outlier.size = 0.5) + 
    theme(
      text = element_text(family = "Arial", size = unit(8, "pt")),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, family = "Arial", size = unit(10, "pt"), face = "bold"),
      legend.text = element_text(family = "Arial", size = unit(8, "pt")), 
      legend.key.size = unit(4, 'pt'),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      axis.text = element_text(family = "Arial", size = unit(8, "pt")),
      aspect.ratio = 1)})
plot <- wrap_plots(plot, ncol = 5)
ggsave(plot, filename = "Plots/Vlnplot-QC.pdf", width = 12, height = 3)

VlnPlot(MergeObject, features = c("nCount_peaks"), group.by = "harmony_clusters", pt.size =0, log = T)

# Annotation
MergeObject$predicted.id = metadata[Cells(MergeObject),]$predicted.id
p1 <- DimPlot(MergeObject, group.by = "Sample", pt.size = 0.1) + DimPlot_theme
p2 <- DimPlot(MergeObject, group.by = "predicted.id", pt.size = 0.1) + DimPlot_theme
p3 <- DimPlot(MergeObject, group.by = "harmony_clusters", pt.size = 0.1, cols = cols) + DimPlot_theme
ggsave(p1 + p2 + p3, filename = "Plots/2.UMAP-sample-clusters.svg", width = 14, height = 4)

p4 <- SpatialPlotV(MergeObject, group.by = "harmony_clusters", ncol = 1, pt.size.factor = 1.1,
  cols = cols, subtitle = FALSE, theme = Spatial_theme + theme(legend.position = "none"))
ggsave(p4, filename = "Plots/2.spatial-cluster.png", width = 5, height = 15, bg = "transparent")
ggsave(p4, filename = "Plots/2.spatial-cluster.svg", width = 5, height = 15)

p5 <- SpatialPlotV(MergeObject, group.by = "predicted.id", ncol = 1,  subtitle = FALSE, 
    theme = Spatial_theme + theme(legend.position = "none"))
ggsave(p5, filename = "Plots/2.spatial-predicted.id.png", width = 2, height = 10, bg = "transparent")
ggsave(p5, filename = "Plots/2.spatial-predicted.id.svg", width = 2, height = 10, bg = "transparent")


plots = list()
for (feature in c("nFrags", "nCount_peaks", "FRIP", "TSSEnrichment", "BlacklistRatio")) {
plots[[feature]] <- SpatialPlotV(MergeObject, feature = feature, subtitle = FALSE, 
    theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1)))
}
combined_plot <- wrap_plots(plots, nrow = 1)
ggsave(combined_plot, filename = "Plots/Spatial-QC.svg", width = 4, height = 4, bg = "transparent")


## Annotations
df = as.data.frame(table(metadata$predicted.id))
types <- df[df$Freq >= 20, ]$Var1
plots.list = c()
for(i in types){
  colors <- rep("#c6c4c4", length(unique(metadata$predicted.id)))
  names(colors) <- unique(metadata$predicted.id)
  colors[i] <- "green" 
  plots.list[[i]] <- SpatialPlotX(MergeObject, group.by = "predicted.id", ncol = 1, cols = colors, subtitle = FALSE,
   theme = Spatial_theme + theme(legend.position = "none")) +
    plot_annotation(title = i)
}
ggsave(cowplot::plot_grid(plotlist = plots.list, nrow = 2), filename = "Plots/annotation.svg", width = 40, height = 25)

df2 = read.csv("Save/stereo.seq.annotations.csv", row.names = 1)
MergeObject$stereo_anno <- df2[Cells(MergeObject),]$annotation
count = data.frame(table(MergeObject$stereo_anno))
plots <- list()
for( i in count[count$Freq > 2,]$Var1){
    highlight <- rownames(MergeObject@meta.data[MergeObject@meta.data$stereo_anno == i,])
    plot <- DimPlot(MergeObject, group.by = "stereo_anno", cells.highlight = highlight,sizes.highlight = 0.8)+ ggtitle(i) + NoLegend() 
    plots[[i]] <- plot
}
plots <- lapply(plots, function(p) {
    p + theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 1
    )})
plot <- wrap_plots(plots, ncol = 6)
ggsave(plot, filename = "Plots/2.stereo-seq_annotation.pdf", width = 20, height = 20)

## Spatial Markergenes ans peaks
MergeObject <- readRDS("Save/MergedObject.rds")
theme_set(theme_void() + custom_theme)
# En2 hs1418 chr5:28178923-28180422
En2.vista <- plot_Spatial_peaks(MergeObject, query_region = "chr5:28178923-28180422")
# Emx2 hs1032 chr19:59465505-59466798
Emx2.vista <- plot_Spatial_peaks(MergeObject, query_region = "chr19:59465505-59466798")
# Hoxa1 hs629 chr6:52319082-52320284
Hoxa1.vista <- plot_Spatial_peaks(MergeObject, query_region = "chr6:52319082-52320284")
# Dlx6 hs298 chr6:6861655-6862377
Dlx6.vista <- plot_Spatial_peaks(MergeObject, query_region = "chr6:6861655-6862377")
# Zic4 hs1043 chr9:91366859-91368203
Zic4.vista <- plot_Spatial_peaks(MergeObject, query_region = "chr9:91366859-91368203")

plots = list()
DefaultAssay(MergeObject) <- "genescore"
for (feature in c("En2", "Emx2", "Dlx6", "Zic4")) {
  plots[[feature]] <- SpatialPlotV(MergeObject, features = feature,  ncol = 1,
    subtitle = FALSE)}

combined_plot <- wrap_plots(
          plots[[1]], En2.vista, 
          plots[[2]], Emx2.vista, 
          plots[[3]], Dlx6.vista,
          plots[[4]], Zic4.vista,
          nrow = 1
          )
ggsave(combined_plot, file = "Plots/SpatialPlot_Markergenes_vista.svg", width = 6, height = 4, bg = "transparent")

SpatialPlotV(subset(MergeObject, 
  subset = (harmony_clusters %in% c("C10","C14"))), group.by = "harmony_clusters", ncol = 1,
    crop = F, image = 1, pt.size.factor = 1, subtitle = FALSE)
#####################################

ggsave(wrap_plots(plots, nrow = 1), filename = "Plots/MarkerGS.1-11.svg", width = 6, height = 4)

labelMarkers = c(
    "Slc4a1", # erythrocyte development
    "Hoxb6", # embryonic skeletal system development
    "Foxc2",  # mesenchyme morphogenesis
    "Nkx2-5",  # heart looping
    "Fgf21", 
    "Fendrr", 
    "Trpv1", 
    "Prrx2",
    "Foxd1"
    )
plots = list()
for (feature in labelMarkers) {
  plots[[feature]] <- SpatialPlotV(MergeObject, features = feature, ncol = 1,
    crop = F, image = 0, pt.size.factor = 1, subtitle = FALSE,
    theme = Spatial_theme + theme(legend.position = "top",
        legend.text = element_text(angle = 45, hjust = 1)))}
Tissues <- wrap_plots(plots, nrow = 1)
ggsave(Tissues, filename = "Plots/MarkerGS.12-20.svg", width = 6, height = 4)

## Genome Tracks
suppressMessages(library(ArchR))
outputDirectory = "Save/MergedObject"
ArchRProj <- loadArchRProject(outputDirectory, showLogo = F)
vista.bed <- rtracklayer::import("/home/disk/zhuming/projects/spatial_ct/public/VISTA/vista.bed")
plot_gene_track <- function(gene, upstream, downstream) {
  p <- plotBrowserTrack(
    ArchRProj = ArchRProj,
    geneSymbol = gene,
    features = list(getPeakSet(ArchRProj), vista.bed),
    loops = getPeak2GeneLinks(ArchRProj),
    groupBy = "harmony_clusters",
    useGroups = paste0("C", 1:11),
    sizes = c(10, 1, 1, 2),
    pal = cols,
    upstream = upstream,
    downstream = downstream,
    normMethod = "ReadsInPromoter",
    threads = getArchRThreads(),
    title = gene,
    borderWidth = 0
  )
  return(p)
}
p2g <- getPeak2GeneLinks(
    ArchRProj = ArchRProj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = F
)
p <- plotPeak2GeneHeatmap(ArchRProj, groupBy = "harmony_clusters")
p
plotPDF(plotList = plot_gene_track("En2", upstream = 120000, downstream = 100000), 
    name = "Plot-Tracks-Marker-En2.pdf", 
    addDOC = FALSE, width = 5, height = 5)

plotPDF(plotList = plot_gene_track("Emx2", upstream = 70000, downstream = 100000), 
    name = "Plot-Tracks-Marker-Emx2.pdf", 
    addDOC = FALSE, width = 5, height = 5)

plotPDF(plotList = plot_gene_track("Zic4", upstream = 30000, downstream = 50000), 
    name = "Plot-Tracks-Marker-Zic4.pdf", 
    addDOC = FALSE, width = 5, height = 5)

plotPDF(plotList = plot_gene_track("Dlx6", upstream = 30000, downstream = 30000), 
    name = "Plot-Tracks-Marker-Dlx6.pdf", 
    addDOC = FALSE, width = 5, height = 5)

plotPDF(plotList = plot_gene_track("Hoxb8", upstream = 100000, downstream = 100000), 
    name = "Plot-Tracks-Marker-Hoxb8.pdf", 
    addDOC = FALSE, width = 4, height = 4)


## Forebrain
DORC.genes <- c("Emx1", "Neurod2", "Neurog2", "Otx1", "Fezf2", "Wnt7b", "Neurod1", "Eomes", "Pax6")
theme_set(theme_void()+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
plots = list()
DefaultAssay(MergeObject) <- "genescore"
for (feature in DORC.genes) {
  plots[[feature]] <- SpatialPlotV(MergeObject, features = feature, subtitle = FALSE,
    theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1))
      )}
combined_plot <- wrap_plots(plots, nrow = 1)
ggsave(combined_plot, file = "Plots/SpatialPlot_Markergenes_Forebrain.svg", width = 9, height = 4, bg = "transparent")

plots = list()
DefaultAssay(MergeObject) <- "RNA"
for (feature in DORC.genes) {
  plots[[feature]] <- SpatialPlotV(MergeObject, features = feature, subtitle = FALSE,
    theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1)),
      cols = c("#384F7F", "#f1e9d8", "#ae2f0f"))}
combined_plot <- wrap_plots(plots, nrow = 1)
ggsave(combined_plot, file = "Plots/SpatialPlot_Genes_Forebrain.svg", width =9, height = 4, bg = "transparent")


## Forebrain Spatial-Temporal Tracks
ArchRProj$Clusters2 <- paste0(ArchRProj$Sample, "_", ArchRProj$harmony_clusters)
colors = c("#EB746A", "#9EA020", "#2CB178", "#2BA5DF", "#B274AF")
names(colors) = paste0(unique(ArchRProj$Sample), "_C4")
peaks  <- getPeakSet(ArchRProj)
# reduced.peaks <- peaks %>%
#   resize(width(peaks) + 1000, fix = "center") %>%
#   reduce(min.gapwidth = 500)
ATAC.peaks <- rtracklayer::import("/home/disk/zhuming/projects/spatial_ct/public/cistromdb/ENCFF838UBQ.bed",format='narrowPeak') # E14.5 Forebrain ATAC
p1 <- plotBrowserTrack(
    ArchRProj = ArchRProj,
    geneSymbol = "Neurod2",
    groupBy = "harmony_clusters",
    pal = cols, 
    sizes = c(15, 1, 1, 1),
    upstream = 70000,
    downstream = 50000,
    normMethod = "ReadsInPromoter",
    threads = getArchRThreads(),
    borderWidth = 0.1,
    tickWidth = 0.1,
    baseSize = 5,
    facetbaseSize = 5,
    plotSummary = c("bulkTrack", "featureTrack"),
    title = "Neurod2"
)
p2 <- plotBrowserTrack(
    ArchRProj = ArchRProj,
    geneSymbol = "Neurod2",
    loops = getCoAccessibility(ArchRProj, resolution = 2000, corCutOff = 0.6),
    groupBy = "Clusters2",
    useGroups = paste0(unique(ArchRProj$Sample), "_C4"),
    features = ATAC.peaks,
    pal = colors, 
    sizes = c(5, 1, 2, 2),
    upstream = 70000,
    downstream = 50000,
    normMethod = "ReadsInPromoter",
    threads = getArchRThreads(),
    borderWidth = 0,
    tickWidth = 0,
    baseSize = 5,
    facetbaseSize = 5,
    plotSummary = c("bulkTrack", "featureTrack"),
    title = "Neurod2"
)
p3 <- plotBrowserTrack(
    ArchRProj = ArchRProj,
    geneSymbol = "Neurod2",
    loops = getCoAccessibility(ArchRProj, resolution = 2000, corCutOff = 0.6),
    groupBy = "Clusters2",
    useGroups = paste0(unique(ArchRProj$Sample), "_C4"),
    pal = colors, 
    sizes = c(1, 1),
    upstream = 70000,
    downstream = 50000,
    normMethod = "ReadsInPromoter",
    threads = getArchRThreads(),
    borderWidth = 0,
    tickWidth = 0,
    baseSize = 5,
    facetbaseSize = 5,
    plotSummary = c("loopTrack", "geneTrack"),
    title = "Neurod2"
)
plotPDF(plotList = p1,
    name = "Plot-Tracks-Marker-Neurod2-Clusters.pdf",
    addDOC = FALSE, width = 4, height = 4)

plotPDF(plotList = p2,
    name = "Plot-Tracks-Marker-Neurod2-Sample.pdf",
    addDOC = FALSE, width = 4, height = 2)

plotPDF(plotList = p3,
    name = "Plot-Tracks-Marker-Neurod2-Loops.pdf",
    addDOC = FALSE, width = 5, height = 1)

##
combine_plots <- plot_Spatial_peaks(MergeObject, "chr11:98259644-98379645", 
  intersect = T, intersect.peaks = reduced.peaks)
ggsave(combine_plots, filename = "Plots/Spatial_Temporal_enhancers.svg", width = 30, height = 8)

## motifs
matches <- getMatches(ArchRProj = ArchRProj, name = "Motif")


plots = list()
DefaultAssay(MergeObject) <- "motif"

abc_genes <- grep("^Tbr1|^Tbr|^Pax6|Neurod2", rownames(MergeObject), value = TRUE)
for (feature in c("Pax6-614", "Tbr1-769", "Neurod2-69")) {
  plots[[length(plots)+1]] <- SpatialPlotV(MergeObject, features = feature,  ncol = 1,
  theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1)),
    subtitle = FALSE, col = rev(RColorBrewer::brewer.pal(9, "RdBu")))}

DefaultAssay(MergeObject) <- "RNA"
for (feature in c("Pax6", "Tbr1", "Neurod2")) {
  plots[[length(plots)+1]] <- SpatialPlotV(MergeObject, features = feature,  ncol = 1,
  theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1)),
    subtitle = FALSE, cols = c("#2968A5", "#EC7085",  "#F4EE63"))}

DefaultAssay(MergeObject) <- "genescore"
for (feature in c("Pax6", "Tbr1", "Neurod2")) {
  plots[[length(plots)+1]] <- SpatialPlotV(MergeObject, features = feature,  ncol = 1,
  theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1)),
    subtitle = FALSE)}

ggsave(wrap_plots(plots, nrow = 1), filename = "Plots/Tbr1.svg", width = 4, height = 4)

