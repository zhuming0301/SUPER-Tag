# rGREAT
setwd("projects/spatial_ct/fig/fig3")
library(rGREAT)
library(clusterProfiler)
suppressMessages(library(dplyr))
suppressMessages(extrafont::loadfonts())

markerPeaks <- readRDS("Save/markerpeaks.rds")
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)
topMarkers <- lapply(markerList, function(df) {
    df <- df[order(-df$Log2FC), ]
    df <- head(df, 20)})
results = list()
get_results_parallel <- function(i) {
  res <- great(topMarkers[[i]], "GO:BP", "txdb:mm10")
  return(getEnrichmentTable(res))
}
results <- mclapply(1:length(topMarkers), get_results_parallel, mc.cores = detectCores())

result <- data.frame()
for (i in 1:length(results)) {
    if (nrow(results[[i]]) > 0) { 
        results[[i]]$cluster <- names(topMarkers)[[i]]
        result <- rbind(result, results[[i]])
    }
}
write.csv(result, 'Save/MarkerPeaks_rGreat_GO_BP_20.csv')

result <- result[order(result$p_adjust_hyper, -result$fold_enrichment_hyper), ]
as.data.frame(result)%>%
    group_by(cluster) %>%
    slice_head(n = 3) %>%
    ungroup() -> top5
top5$logp = -log10(top5$p_adjust_hyper)
top5$cluster = factor(top5$cluster, levels = paste0("C",sort(as.numeric(gsub("C", "", unique(top5$cluster))))))
top5 <- top5[order(top5$cluster, top5$logp), ]
top5$description = factor(top5$description, levels = rev(unique(top5$description)))
p <- ggplot(top5, aes(x=cluster,y=description, color=logp, size = fold_enrichment)) + 
    geom_point(stat="identity") +
    ggtitle("GO enrichment of markerPeaks") + 
    scale_size(range = c(1, 3)) +
    guides(
    color = guide_colorbar(title = "-log10(p adjust hyper)", 
                           barwidth = 1, 
                           barheight = 5, 
                           override.aes = list(size = 2)),
    size = guide_legend(title = "fold enrichment")
  ) + ylab("") + xlab("")
ggsave(p, filename = "Plots/GO.MarkerPeaks.top20.pdf", width = 4.5, height = 3)





