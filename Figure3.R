library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readxl)
library(ggdendro)
library(patchwork)

filep <- "./MiddleFile/"

group_order <- c("Control", "Alive", "Deceased")
group_colors <- c(
  Control  = "#7ac7e2",
  Alive    = "#f9d580",
  Deceased = "#e3716e"
)

module_order <- c(
  "IL6_JAK_STAT3",
  "TNF_NFkB",
  "Inflammatory_response",
  "IFN_gamma",
  "Complement",
  "Coagulation",
  "Antigen_presentation",
  "Cytotoxicity",
  "ROS_pathway",
  "Glycolysis",
  "OXPHOS"
)

ct_order_main <- c("Mono", "DC", "B", "CD4+T", "CD8+T", "NK")
ct_order_augur <- c("Plasma", "Mono", "DC", "CD8+T", "B", "NK", "CD4+T")

pub_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold")
  )

pb <- readRDS(paste0(filep, "Fig3_pseudobulk_revised.rds"))
donor_meta <- pb$donor_meta

## Figure 3A
pca_df <- read_xlsx(
  paste0(filep, "Fig3A_PCA_SourceData.xlsx"),
  sheet = "PCA_plot_data"
) %>%
  mutate(
    outcome = factor(outcome, levels = group_order)
  )

p3A_variance <- read_xlsx(
  paste0(filep, "Fig3A_PCA_SourceData.xlsx"),
  sheet = "PCA_variance"
)
percent_var <- p3A_variance$Percent_variance

p3A_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = outcome, fill = outcome)
) +
  geom_point(
    size = 3.2,
    alpha = 0.9,
    shape = 21,
    stroke = 0.35,
    color = "black"
  ) +
  stat_ellipse(
    aes(color = outcome),
    type = "norm",
    level = 0.68,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_blank()
  )

p3A_pca

## Figure 3B
pca_centroids <- pca_df %>%
  group_by(outcome) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("outcome")

dist_conditions <- dist(pca_centroids)
hc <- hclust(dist_conditions)

hc_dendro <- dendro_data(hc)

p3_hc <- ggdendrogram(hc_dendro) +
  geom_text(
    data = hc_dendro$labels,
    aes(x = x, y = y, label = label, color = label),
    vjust = 1.5,
    fontface = "bold",
    size = 4
  ) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.2, 0)) +
  ylab("Height") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "none"
  )

p3_hc

## Figure 3C
ss_donor <- readRDS(paste0(filep, "Fig3_ssGSEA_donor_revised.rds"))
kw_3b <- read.csv(paste0(filep, "Figure3B_pathway_KW_BH_stats_11pathways.csv"))

dl <- as.data.frame(t(ss_donor)) %>%
  rownames_to_column("orig.ident") %>%
  left_join(donor_meta, by = "orig.ident") %>%
  pivot_longer(
    cols = all_of(rownames(ss_donor)),
    names_to = "module",
    values_to = "score"
  ) %>%
  mutate(
    outcome = factor(outcome, levels = group_order),
    module = factor(module, levels = module_order)
  )

gm <- dl %>%
  group_by(module, outcome) %>%
  summarise(
    grp_mean = mean(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(module) %>%
  mutate(
    z = as.numeric(scale(grp_mean))
  ) %>%
  ungroup() %>%
  mutate(
    module = factor(module, levels = rev(module_order)),
    outcome = factor(outcome, levels = group_order)
  )

lab_df <- kw_3b %>%
  mutate(
    label = ifelse(
      sig == "ns",
      as.character(module),
      paste0(module, "  ", sig)
    )
  )

lab_map <- setNames(lab_df$label, lab_df$module)

p3C <- ggplot(gm, aes(x = outcome, y = module, fill = z)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(
    low = "#7ac7e2",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    name = "Z-score"
  ) +
  scale_y_discrete(labels = lab_map) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_fixed() +
  pub_theme +
  theme(
    plot.title = element_blank(),
    plot.caption = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

p3C

## Figure 3D
ct_pathway <- read_xlsx(
  paste0(filep, "Fig3D_Pathway_Celltype_SourceData.xlsx"),
  sheet = "Pathway_celltype_summary"
) %>%
  mutate(
    module = factor(module, levels = rev(module_order)),
    celltype_major = factor(celltype_major, levels = ct_order_main)
  )

p3D <- ggplot(
  ct_pathway,
  aes(x = celltype_major, y = module, fill = z)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(
    low = "#7ac7e2",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    name = "Z-score"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_fixed() +
  pub_theme +
  theme(
    plot.title = element_blank(),
    plot.caption = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

p3D

## Figure 3E
augur_Control_vs_Alive <- readRDS(paste0(filep, "Fig3A_augur_celltype_major_Control_vs_Alive.rds"))
augur_Control_vs_Deceased <- readRDS(paste0(filep, "Fig3A_augur_celltype_major_Control_vs_Deceased.rds"))
augur_Alive_vs_Deceased <- readRDS(paste0(filep, "Fig3A_augur_celltype_major_Alive_vs_Deceased.rds"))

plot_augur_auc <- function(
    augur_res,
    title = "",
    colors = c("#7ac7e2", "#fabb6e"),
    ct_order = ct_order_augur
) {

  data <- augur_res[["AUC"]] %>%
    as.data.frame() %>%
    mutate(
      cell_type = factor(cell_type, levels = rev(ct_order))
    )

  ggplot(data, aes(x = auc, y = cell_type, color = auc)) +
    geom_segment(
      aes(x = 0, xend = auc, y = cell_type, yend = cell_type),
      linewidth = 0.5,
      color = "grey55"
    ) +
    geom_point(size = 3) +
    geom_vline(
      xintercept = mean(data$auc, na.rm = TRUE),
      linetype = "dotted",
      color = "black",
      linewidth = 0.5
    ) +
    scale_color_gradient(
      low = colors[1],
      high = colors[2],
      name = "AUC"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      title = title,
      x = "AUC",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

pa1 <- plot_augur_auc(
  augur_Control_vs_Alive,
  title = "Control vs Alive",
  colors = c("#7ac7e2", "#fad354")
)

pa2 <- plot_augur_auc(
  augur_Control_vs_Deceased,
  title = "Control vs Deceased",
  colors = c("#7ac7e2", "#b7282e")
)

pa3 <- plot_augur_auc(
  augur_Alive_vs_Deceased,
  title = "Alive vs Deceased",
  colors = c("#fad354", "#b7282e")
)

p_augur <- pa1 | pa2 | pa3
p_augur
