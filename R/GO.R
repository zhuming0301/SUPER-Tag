library(org.Mm.eg.db)
library(GOSemSim)
library(parallel)
library(ggplot2)
library(clusterProfiler)
library(ArchR)
library(dplyr)
suppressMessages(extrafont::loadfonts())
custom_theme <- theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 6),  # 全局字体设置
    plot.title = element_text(family = "Arial", size = 7, face = "bold"),  # 标题字体
    axis.title = element_text(family = "Arial", size = 6),  # 坐标轴标题字体
    axis.text = element_text(family = "Arial", size = 5),  # 坐标轴刻度字体
    legend.text = element_text(family = "Arial", size = 5),  # 图例字体
    legend.key.size = unit(1, "pt")
  )
theme_set(custom_theme)
markerGS <- readRDS(file = "Save/MarkerGeneScore.11-20.rds")
markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
topMarkers <- lapply(markerList, function(df) {
    df <- df[order(-df$Log2FC), ]
    head(df, 50)})
geneLists <- lapply(topMarkers, function(x) x$name) # 提取所有分组的基因名列表
performGO <- function(genes) {
  enrichGO(
    gene = genes, 
    OrgDb = org.Mm.eg.db, 
    keyType = "SYMBOL", 
    ont = "BP", # 可选 "BP", "MF", "CC"
    pAdjustMethod = "BH", 
    qvalueCutoff = 0.05
  )
}
# 并行执行
goResults <- mclapply(geneLists, performGO, mc.cores = detectCores()/2)
names(goResults) <- names(markerList)

validResultsWithNames <- lapply(names(goResults), function(clusterName) {
  # 获取当前分组的 GO 富集结果
  result <- goResults[[clusterName]]
  # 如果是有效数据，添加一列 Cluster 并返回数据框
  if (!is.null(result) && nrow(result) > 0) {
    result_df <- as.data.frame(result)  # 转换为数据框
    result_df$Cluster <- clusterName  # 添加分组名称列
    return(result_df)  # 返回带有分组名称的结果
  } 
})

# 过滤掉 NULL 元素
validResultsWithNames <- validResultsWithNames[!sapply(validResultsWithNames, is.null)]
mergedResults <- do.call(rbind, validResultsWithNames)
write.csv(mergedResults, file = "Save/MarkerGeneScore_GO_BP.12-20.csv", row.names = FALSE)

# 准备数据矩阵
mergedResults <- read.csv("Save/MarkerGeneScore_GO_BP.1-11.csv")
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
                           override.aes = list(size = 2)),  # 调整 color 图例的点大小
    size = guide_legend(title = "Gene Ratio")
  ) + scale_size(range = c(0.5, 3)) + custom_theme
ggsave(p, filename = "Plots/GO.MarkerGS.top.12-20.pdf", width = 4.5, height = 3)




