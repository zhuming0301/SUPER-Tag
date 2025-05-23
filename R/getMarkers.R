setwd("projects/spatial_ct/")
set.seed(8)
suppressMessages(library(ArchR))
addArchRGenome("mm10")
addArchRThreads(threads = 24)
source("Code/ArchRtoSignac.R")
outputDirectory = "Save/MergedObject"
ArchRProj <- loadArchRProject(outputDirectory)

# Marker peaks
markerPeaks <- getMarkerFeatures(
    ArchRProj, 
    useGroups = paste0("C", 1:11),
    bgdGroups = paste0("C", 12:20),
    useMatrix = "PeakMatrix", 
    groupBy = "harmony_clusters"
    )
heatmapPeaks <- plotMarkerHeatmap(
    markerPeaks, 
    cutOff = "FDR <= 0.01 & Log2FC >= 1", 
    limits = c(-3, 3), nLabel = 1, clusterCols = F,
    )
pdf("Plots/MarkerPeaks2.pdf", width = 6, height = 8) 
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveRDS(markerPeaks, file = "Save/markerpeaks.rds")

# Maker genescore
markerGS <- getMarkerFeatures(
    ArchRProj, 
    useMatrix = "GeneScoreMatrix", groupBy = "harmony_clusters", 
    useGroups = c(paste0("C", 1:11)),
    bgdGroups = c(paste0("C", 12:20))
    )

labelMarkers = c(
    "Dbx1", "En1", "En2", # Mb Vz
    "Lhx9", "Tfap2b", "Pou4f1", # Mb
    "Lhx2", "Emx1", # Pall VZ
    "Neurod6", "Tbr1", # Pall
    "Dlx5", "Dlx1", # SPall 
    "Igfbpl1", "Lhx5","Otx1","Acsl6", 
    "Pax6","Hoxb8","Neurog3", "Neurod4","Neurod2","Ina", "Zic4","Pax3",
    "Hoxc5", "Prdm13", "Hoxa1", "Pgpep1", "Foxg1", "Emx2",
    "Six3", "Olig2", "Olig1", "Dlx2", "Gad2", "Gad1", "Dlx6",
    "Myo3b", "Lhx1", "Myo18b", "Kcnq2", "Nr2f1", 
    "Tfap2d", "Lhx9", "Pou4f2", "Olfr133", "Vmn2r49", "Olfr123"
    )

markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

gene_names <- rowData(markerGS)$name
keep <- !grepl("Rik|^Gm|^mt-|^Rps|^Rpl|^Mir", gene_names)
markerGS_filtered <- markerGS[keep, ]

topMarkers <- lapply(markerGS_filtered, function(df) {
    df <- df[order(-df$Log2FC), ]
    head(df, 50)})
selectedPeaks <- unique(unlist(lapply(topMarkers, function(df) rownames(df))))
seMarker_subset <- markerGS_filtered[selectedPeaks, ]
heatmapGS <- plotMarkerHeatmap(seMarker_subset, limits = c(-3, 3), nLabel = 1, clusterCols = F)

heatmapGS <- plotMarkerHeatmap(
    markerGS_filtered, 
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
    limits = c(-3, 3), nLabel = 10, clusterCols = F,
    labelMarkers = labelMarkers
    )

pdf("Plots/MarkerGeneScore.pdf", width = 5, height = 7) 
draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveRDS(markerGS_filtered, file = "Save/MarkerGeneScore.rds")
# markerGS_filtered <- readRDS("Save/MarkerGeneScore.rds")


# Maker genescore tissue

markerPeaks2 <- getMarkerFeatures(
    ArchRProj, 
    bgdGroups = paste0("C", 1:11),
    useGroups = paste0("C", 12:20),
    useMatrix = "PeakMatrix", 
    groupBy = "harmony_clusters"
    )
heatmapPeaks <- plotMarkerHeatmap(
    markerPeaks2, 
    cutOff = "FDR <= 0.01 & Log2FC >= 1", 
    limits = c(-3, 3), nLabel = 1, clusterCols = F,
    )
pdf("Plots/MarkerPeaks.12-20.pdf", width = 6, height = 8) 
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveRDS(markerPeaks2, file = "Save/markerpeaks2.rds")


markerGS2 <- getMarkerFeatures(
    ArchRProj, 
    useMatrix = "GeneScoreMatrix", groupBy = "harmony_clusters", 
    bgdGroups = c(paste0("C", 1:11)),
    useGroups = c(paste0("C", 12:20))
    )
gene_names <- rowData(markerGS2)$name
keep <- !grepl("Rik|^Gm|^mt-|^Rps|^Rpl|^Mir", gene_names)
markerGS_filtered2 <- markerGS2[keep, ]

labelMarkers = c(
    "Klf1", "Slc4a1", "Hba-a2", # erythrocyte development
    "Six1", "Hoxb6", "Hoxc5", # embryonic skeletal system development
    "Hoxb8", "Ccl7", 
    "Foxf1", "Foxc2",  # mesenchyme morphogenesis
    "Hand2", "Nkx2-5", "Gata4", # heart looping
    "Fgf21", "Hcar1", "Tlx2", # Mb
    "Fendrr", "Gata5", 
    "Atp2a1", "Trpv1", "Tnni2",
    "Prrx2", "Msx1", "Osr1", "Alx4",
    "Matn4", "Foxd1"
    )

markerList2 <- getMarkers(markerGS_filtered2, cutOff = "FDR <= 0.01 & Log2FC >= 1")
topMarkers <- lapply(markerList2, function(df) {
    df <- df[order(-df$Log2FC), ]
    head(df, 50)})
selectedPeaks <- unique(unlist(lapply(topMarkers, function(df) rownames(df))))
seMarker_subset <- markerGS_filtered2[selectedPeaks, ]

heatmapGS <- plotMarkerHeatmap(
    seMarker_subset, 
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
    limits = c(-3, 3), nLabel = 5, clusterCols = F,
    labelMarkers = labelMarkers
    )

pdf("Plots/MarkerGeneScore.12-20.pdf", width = 4, height = 6) 
draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

saveRDS(markerGS_filtered2, file = "Save/MarkerGeneScore.11-20.rds")
# markerGS_filtered2 <- readRDS("Save/MarkerGeneScore.11-20.rds")



## Motif
ArchRProj <- addMotifAnnotations(ArchRProj = ArchRProj, motifSet = "cisbp", name = "Motif")
markerPeaks <- readRDS("Save/markerpeaks.rds")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = ArchRProj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs, rastr = F, n = 15)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 6, height = 8, addDOC = FALSE)
## Plotting ComplexHeatmap!
ArchRProj <- addBgdPeaks(ArchRProj)
ArchRProj <- addDeviationsMatrix(
  ArchRProj = ArchRProj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(ArchRProj, name = "MotifMatrix", plot = F)

motifs <- c("Neurod1", "Sox2", "Eomes", "Tbr1", "Pax6")
markerMotifs <- getFeatures(ArchRProj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

p <- plotEmbedding(
    ArchRProj = ArchRProj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(ArchRProj)
)
p2 <- lapply(p, function(x){
    x + guides(color = "none", fill = "none") + 
    theme_void() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

## get motifs
# Neurod2 "chr11:98259644-98379645"
pSet <- getPeakSet(ArchRProj = ArchRProj)
pSet <- pSet %>%
  resize(width(pSet) + 1000, fix = "center") %>%
  reduce(min.gapwidth = 500)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")

matches <- getMatches(ArchRProj = ArchRProj, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(98259644), end = c(98379645)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
tf.list = c()
for(i in queryHits){
    tf.list[[length(tf.list) + 1]] = colnames(matches)[which(assay(matches[i,]))]
}

seGroupMotif <- getGroupSE(ArchRProj, useMatrix = "MotifMatrix", groupBy = "harmony_clusters")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGIM_MM <- correlateMatrices(ArchRProj,
    k=10,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)


corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
write.csv(corGIM_MM, file = "Save/Motifs.csv", quote = F, row.names = F)
TFRegulators <- sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
tf.list <- lapply(tf.list, function(x){x = intersect(sub("_.*", "",x), TFRegulators)})
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p



## plot Spatial motif
MergeObject <- readRDS("Save/MergedObject.genescore.rds")
MergeObject <- addMotifAssay(ArchRProj, MergeObject)
DefaultAssay(MergeObject) <- "motif"

motifs <- c("Neurod1", "Sox2", "Eomes", "Tbr1", "Pax6")

total.motifs <- rownames(MergeObject[["motif"]])
matched <- unlist(sapply(motifs, function(x) grep(paste0("^", x), total.motifs, value = TRUE)))

plots = list()
for (feature in matched) {
  plots[[feature]] <- SpatialPlotV(MergeObject, features = feature,  ncol = 1,
    subtitle = FALSE, theme = Spatial_theme + theme(
      legend.position = "top",
      legend.text = element_text(angle = 45, hjust = 1)),
      cols = ArchRPalettes$solarExtra)}

combined_plot <- wrap_plots(plots, nrow = 1)
ggsave(combined_plot, file = "Plots/Spatial-Motif.png", width = 30, height = 5)
