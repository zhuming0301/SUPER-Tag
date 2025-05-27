# 
SpatialPlotH <- function(
    object, 
    group.by = NULL, 
    features = NULL,
    pt.size.factor = 1, 
    crop = FALSE, 
    image.alpha = 0, 
    cols = NULL, 
    max.cutoff = 'q99',
    subtitle = TRUE,
    legend = TRUE,
    legend_rel_widths = c(15, 1),
    legend_ncol = 2,
    plot_ncol = 2,
    ncol = NULL,
    ...
    ) {
        plots <- SpatialPlot(
            object = object, 
            group.by = group.by, 
            features = features, 
            pt.size.factor = pt.size.factor, 
            cols = cols,
            crop = crop,  
            image.alpha = image.alpha, 
            stroke = 0, 
            max.cutoff = max.cutoff, 
            combine = FALSE, 
            ...
            )
        if(subtitle){plots_no_legend <- lapply(plots, function(p) {
            # p$layers[[1]]$aes_params$shape=22
            p + theme(legend.position = "none")}
            )} else {
               plots_no_legend <- lapply(plots, function(p) {
            # p$layers[[1]]$aes_params$shape=22
            p + theme(legend.position = "none") + ggtitle("")}
            )
            }

        if(is.null(ncol)){
          ncol = length(plots_no_legend)
        }
        if(legend){
          plots_no_legend[[length(plots_no_legend)]] <- plots_no_legend[[length(plots_no_legend)]] + theme(legend.position = "right")
          final_plot <- wrap_plots(plots_no_legend, ncol = ncol)
        }
         else {
          final_plot <- wrap_plots(plots_no_legend, ncol = ncol)
        }
    final_plot
}

# 
SpatialPlotV <- function(
    object, 
    group.by = NULL, 
    features = NULL,
    pt.size.factor = 1, 
    crop = FALSE, 
    image.alpha = 0, 
    cols = NULL, 
    max.cutoff = 'q99',
    subtitle = TRUE,
    theme = NULL,
    ...
    ) {
        theme_set(theme_void())
        plots <- SpatialPlot(
            object = object, 
            group.by = group.by, 
            features = features, 
            pt.size.factor = pt.size.factor, 
            cols = cols,
            crop = crop,  
            image.alpha = image.alpha, 
            stroke = 0, 
            max.cutoff = max.cutoff, 
            combine = FALSE, 
            ...
            )
        if(subtitle == FALSE){
            plots <- lapply(plots, function(p) {p + ggtitle("")})
            } 
        if(!is.null(features) && !is.null(cols)){
          plots <- lapply(plots, function(p) {
            p + scale_fill_gradientn(colors = cols)
            })
            }
        plots <- lapply(plots, function(p) {p + theme(legend.position = "none")})
        if(!is.null(theme)){
            plots[[1]] <- plots[[1]] + theme
            } else {
            plots[[1]] <- plots[[1]] + theme(legend.position = "top")
            }
        final_plot <- wrap_plots(plots, ncol = 1)
        final_plot
}

plot_Spatial_peaks <- function(project, query_region, assay = "peaks", intersect = FALSE, intersect.peaks = NULL){
  strs <- strsplit(query_region,  split = ":|-")[[1]]
  query_region <- GenomicRanges::GRanges(seqnames = strs[[1]], ranges = IRanges::IRanges(start = as.integer(strs[[2]]), end = as.integer(strs[[3]])))
  if(intersect){
    if(is.null(intersect.peaks)){
      intersect.peaks <- granges(project)
    }     
    peaks.plot <- GenomicRanges::intersect(query_region, intersect.peaks)
    } else {
    peaks.plot <- query_region
  }
  fm <- FeatureMatrix(
    Fragments(project[[assay]]),
    peaks.plot,
    cells = NULL,
    process_n = 2000,
    sep = c(":", "-"),
    verbose = TRUE)
  a <- t(as.data.frame(fm))
  colnames(a)
  project@meta.data <- cbind(project@meta.data, a)
  plot.list = c()
  for(peak in colnames(a)){
    plot.list[[length(plot.list)+1]] <- SpatialPlotV(project, features = peak, cols = c("#cdcdcd", "red"),
    crop = F, pt.size.factor = 1, subtitle = FALSE, max.cutoff = quantile(fm, probs = 0.98), ncol = 1)
  }
  plot.list
  combine_plots = wrap_plots(plot.list, nrow = 1)
  combine_plots
}
