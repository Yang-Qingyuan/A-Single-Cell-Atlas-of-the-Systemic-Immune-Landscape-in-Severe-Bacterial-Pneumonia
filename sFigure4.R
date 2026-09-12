library(rio)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)

filep <- "./MiddleFile/"

th <- theme(text = element_text(color = "black"),
            axis.title = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 12, face = "bold"),
            plot.caption = element_text(size = 8),
            axis.line = element_line(linewidth = 1),
            axis.ticks = element_line(linewidth = 0.5),
            axis.text = element_text(size = 10, face = "bold", color = "black"))

outl <- c("Control", "Alive", "Deceased")
outcomeCol <- c('#7ac7e2', '#f9d580', '#e3716e')
hm_cols <- colorRampPalette(c('#8ccbea', '#f7df87', '#e68b81'))(100)

## Figure S4A
mono_prop <- import(paste0(filep, "FigS4A.xlsx"))
mono_prop$outcome <- factor(mono_prop$outcome, levels = outl)
mono_prop$celltype <- factor(mono_prop$celltype,
                             levels = c("Mono_c01-Classical", "Mono_c02-Intermediate", "Mono_c03-Non_Classical"))

p4.1 <- ggplot(mono_prop, aes(x = outcome, y = prop, fill = outcome)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    color = "black"
  ) +
  geom_jitter(
    aes(color = outcome),
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(~celltype, scales = "free_y") +
  scale_fill_manual(values = outcomeCol) +
  scale_color_manual(values = outcomeCol) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "",
    y = "Proportion of all cells",
    title = "scRNA-seq"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  )
p4.1

## Figure S4C
s4c <- import_list(paste0(filep, "FigS4C_Monocyte_functional_heatmap_SourceData.xlsx"))

get_mat <- function(df) {
  rownames(df) <- df$Function
  df$Function <- NULL
  as.matrix(df)
}

p_c01 <- pheatmap(
  get_mat(s4c[["C01_heatmap_Zscore"]]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = hm_cols,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Mono c01 Classical",
  angle_col = "90",
  silent = TRUE
)

p_c02 <- pheatmap(
  get_mat(s4c[["C02_heatmap_Zscore"]]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = hm_cols,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Mono c02 Intermediate",
  angle_col = "90",
  silent = TRUE
)

p_c03 <- pheatmap(
  get_mat(s4c[["C03_heatmap_Zscore"]]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = hm_cols,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Mono c03 Non Classical",
  angle_col = "90",
  silent = TRUE
)

grid.arrange(
  p_c01[[4]],
  p_c02[[4]],
  p_c03[[4]],
  ncol = 3
)

## Figure S4F
nk1 <- import_list(paste0(filep, "FigS4F_NK_c01_CD56high_heatmap_SourceData.xlsx"))[["Mean_expression"]]
rownames(nk1) <- nk1$Gene
nk1$Gene <- NULL
hm1 <- as.matrix(nk1)

p4.9 <- pheatmap(hm1,
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 color = hm_cols,
                 angle_col = '90',
                 gaps_row = c(8),
                 scale = 'row',
                 main = 'NK_c01-CD56high')
p4.9

nk2 <- import_list(paste0(filep, "FigS4F_NK_c02_CD56low_heatmap_SourceData.xlsx"))[["Mean_expression"]]
rownames(nk2) <- nk2$Gene
nk2$Gene <- NULL
hm2 <- as.matrix(nk2)

p4.10 <- pheatmap(hm2,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE,
                  color = hm_cols,
                  angle_col = '90',
                  scale = 'row',
                  main = 'NK_c02-CD56low')
p4.10

## Figure S4G
mdc <- import(paste0(filep, "FigS4G_mDC_DotPlot_SourceData.xlsx"))
mdc$Gene <- factor(mdc$Gene, levels = unique(mdc$Gene))
mdc$Outcome <- factor(mdc$Outcome, levels = outl)

p4.8 <- ggplot(mdc, aes(x = Gene, y = Outcome)) +
  geom_point(aes(size = Percent_expressed, color = Scaled_average_expression)) +
  scale_color_gradient(low = '#8ccbea', high = '#e68b81') +
  coord_flip() +
  labs(x = '', y = '', title = 'mDC') +
  theme_classic() +
  th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p4.8
