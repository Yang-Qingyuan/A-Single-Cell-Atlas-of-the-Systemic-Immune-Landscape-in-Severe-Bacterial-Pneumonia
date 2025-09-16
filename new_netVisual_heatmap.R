library(ComplexHeatmap)
library(reshape2)
new_netVisual_heatmap <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                                  title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                                  sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL, color_range = NULL,
                                  max_row_barplot = NULL, max_col_barplot = NULL){
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  
  net <- net.diff
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  
  if (is.null(color_range)) {
    min_val <- min(mat, na.rm = TRUE)
    max_val <- max(mat, na.rm = TRUE)
  } else {
    min_val <- color_range[1]
    max_val <- color_range[2]
  }
  
  if (min_val < 0) {
    color.heatmap.use = colorRamp3(c(min_val, 0, max_val), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min_val, max_val), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min_val, max_val), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
  }
  
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1 <- rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE, gp = gpar(fill = color.use, col=color.use),
                                               ylim = if (!is.null(max_row_barplot)) c(0, max_row_barplot) else NULL), 
                       show_annotation_name = FALSE)
  ha2 <- HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE, gp = gpar(fill = color.use, col=color.use),
                                                   ylim = if (!is.null(max_col_barplot)) c(0, max_col_barplot) else NULL),
                           show_annotation_name = FALSE)
  
  ht1 <- Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                 bottom_annotation = col_annotation, left_annotation = row_annotation, top_annotation = ha2, right_annotation = ha1,
                 cluster_rows = cluster.rows, cluster_columns = cluster.cols,
                 row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), column_names_gp = gpar(fontsize = font.size),
                 column_title = title.name, column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90,
                 row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90,
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                             border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")))
  return(ht1)
}

