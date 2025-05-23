setwd("projects/spatial_ct/")
options(warn = -1)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(ArchR))
suppressMessages(library(parallel))
addArchRThreads(threads = 8)
source("Code/SpatialPlot_new.R")
outputDirectory = "Save/MergedObject"
ArchRProj <- loadArchRProject(outputDirectory, showLogo = F)
metadata <- getCellColData(ArchRProj)
C4.cells <- row.names(metadata[metadata$harmony_clusters == "C4",])
cortex <- subsetArchRProject(ArchRProj, cells = C4.cells, outputDirectory = "Save/C4")

plotVarDev <- getVarDeviations(cortex, name = "MotifMatrix", plot = F)
write.csv(plotVarDev, file = "Save/Motifs.C4.csv", quote = F, row.names = F)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = cortex, addDOC = FALSE)

## Plotting Ggplot!
pSet <- getPeakSet(ArchRProj)
peaks <- pSet %>%
  resize(width(peaks) + 1000, fix = "center") %>%
  reduce(min.gapwidth = 1000)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")

matches <- getMatches(ArchRProj, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]


gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(98259644), end = c(98350000)))
candidates <- GenomicRanges::intersect(gr, peaks)
queryHits <- queryHits(findOverlaps(query = pSet, subject = candidates[3], type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = ArchRProj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 1"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 15)



# 
outputDirectory = "Save/C4"
cortex <- loadArchRProject(outputDirectory, showLogo = F)
ATAC.peaks <- rtracklayer::import("/home/disk/zhuming/projects/spatial_ct/public/cistromdb/ENCFF838UBQ.bed",format='narrowPeak') # E14.5 Forebrain ATAC
gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(98259644), end = c(98350000)))

cortex <- addPeakSet(cortex, intersect(ATAC.peaks, gr), force=TRUE)
cortex <- addMotifAnnotations(ArchRProj = cortex, motifSet = "cisbp", name = "Motif", force = TRUE)

pSet <- getPeakSet(ArchRProj = cortex)
rtracklayer::export(
    pSet, 
    con = "Save/Neurod2_peaks.bed", 
    format = "bed"
    )
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
matches <- getMatches(ArchRProj = cortex, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]


queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))

df.list = list()
for(i in 1:length(pSet)){
  print(i)
  df.list[[i]] = colnames(matches)[which(assay(matches[i,]))]
}
save(df.list, file = "Motifs.rdata")

load("Motifs.rdata")



cortex <- addClusters(cortex, reducedDims = "Harmony")
seGroupMotif <- getGroupSE(cortex, useMatrix = "MotifMatrix", groupBy = "Clusters")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGIM_MM <- correlateMatrices(
    ArchRProj = cortex,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.6 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.6))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES", 1])
label_data <- as.data.frame(corGIM_MM[corGIM_MM$TFRegulator=="YES",])
library(ggrepel) 
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point(size = 1) + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  ) +
  geom_text_repel(
  data = label_data,
  aes(label = GeneIntegrationMatrix_name),
  color = "black",         
  size = 3,                
  box.padding = 0.2,        
  segment.color = "grey50", 
  max.overlaps = Inf
  )
ggsave(p,filename = "Plots/Motifs.cor.pdf", width = 5, height = 5)
