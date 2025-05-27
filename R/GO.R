library(org.Mm.eg.db)
library(GOSemSim)
library(parallel)
library(ggplot2)
library(clusterProfiler)
library(ArchR)
library(dplyr)
suppressMessages(extrafont::loadfonts())

markerGS <- readRDS(file = "Save/MarkerGeneScore.rds")
markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
topMarkers <- lapply(markerList, function(df) {
    df <- df[order(-df$Log2FC), ]
    head(df, 50)})
geneLists <- lapply(topMarkers, function(x) x$name) 
performGO <- function(genes) {
  enrichGO(
    gene = genes, 
    OrgDb = org.Mm.eg.db, 
    keyType = "SYMBOL", 
    ont = "BP",
    pAdjustMethod = "BH", 
    qvalueCutoff = 0.05
  )
}

goResults <- mclapply(geneLists, performGO, mc.cores = detectCores()/2)
names(goResults) <- names(markerList)

validResultsWithNames <- lapply(names(goResults), function(clusterName) {
  result <- goResults[[clusterName]]
  if (!is.null(result) && nrow(result) > 0) {
    result_df <- as.data.frame(result) 
    result_df$Cluster <- clusterName 
    return(result_df) 
  } 
})

validResultsWithNames <- validResultsWithNames[!sapply(validResultsWithNames, is.null)]
mergedResults <- do.call(rbind, validResultsWithNames)
write.csv(mergedResults, file = "Save/MarkerGeneScore_GO_BP.csv", row.names = FALSE)

mergedResults <- read.csv("Save/MarkerGeneScore_GO_BP.csv")
mergedResults <- mergedResults[order(mergedResults$Cluster, mergedResults$p.adjust), ]
as.data.frame(mergedResults)%>%
    group_by(Cluster) %>%
    slice_head(n = 3) %>%
    ungroup() -> top5
top5$logp = -log10(top5$p.adjust)
top5$Cluster = factor(top5$Cluster, levels = paste0("C",sort(as.numeric(gsub("C", "", unique(top5$Cluster))))))
top5 <- top5[order(top5$Cluster, top5$p.adjust), ]
top5$Description = factor(top5$Description, levels = rev(unique(top5$Description)))
top5$GeneRatio = sapply(top5$GeneRatio, function(x) eval(parse(text = x)))
top5$BgRatio = sapply(top5$BgRatio, function(x) eval(parse(text = x)))
p <- ggplot(top5, aes(x=Cluster,y=Description, color=logp, size = GeneRatio)) + geom_point(stat="identity") + theme_bw() +
    theme(axis.text.y = element_text(size = 12)) + 
    ggtitle("GO enrichment of GeneScore markers") + 
    guides(
    color = guide_colorbar(title = "-log10(p adjust)", 
                           barwidth = 1, 
                           barheight = 5, 
                           override.aes = list(size = 2)),
    size = guide_legend(title = "Gene Ratio")
  ) + scale_size(range = c(0.5, 3))
ggsave(p, filename = "Plots/GO.MarkerGS.top.pdf", width = 4.5, height = 3)




