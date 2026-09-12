library(rio)
library(dplyr)
library(ggplot2)
library(ggpubr)

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

## Figure S5A
b_prop <- import(paste0(filep, "FigS5A_Bcell_proportion_allcells_SourceData.xlsx"))
b_prop$outcome <- factor(b_prop$outcome, levels = outl)
b_prop$celltype <- factor(b_prop$celltype,
                          levels = c('B_c01-mature naive B', 'B_c02-memory B',
                                     'B_c03-active Plasma', 'B_c04-terminal Plasma'))

p4.1 <- ggplot(b_prop, aes(x = outcome, y = prop, fill = outcome)) +
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
  facet_wrap(~celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = outcomeCol) +
  scale_color_manual(values = outcomeCol) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "",
    y = "Proportion of all cells",
    title = "scRNA-seq"
  ) +
  theme_classic2() +
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

## Figure S5D
meta_agg <- import(paste0(filep, 'FigS5D.xlsx'))
bcricol <- c("#b5dcf3", "#71b8ec",
             "#c1e2db", "#82c4b8",
             "#f9dbb5", "#f1b56d", '#f9bdb5', '#f67c6f')
meta_agg$BCR_isotype <- factor(meta_agg$BCR_isotype, levels =
                                 c('IGHM', 'IGHD', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'))
meta_agg$outcome <- factor(meta_agg$outcome, levels = outl)

p2 <- ggplot(data = meta_agg, mapping = aes(x = celltype, y = percentage, fill = BCR_isotype)) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.7) +
  scale_fill_manual(values = bcricol) +
  labs(y = "BCR isotype (%)", x = "", title = "", fill = "") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) + th
p2

## Figure S5E
meta_agg <- import(paste0(filep, 'FigS5E.xlsx'))
bcrcol1 <- c("#cbe3e6", "#efc661")
meta_agg$clonestatus <- factor(meta_agg$clonestatus, levels =
                                 c("No clonal", "Clonal"))
meta_agg$outcome <- factor(meta_agg$outcome, levels = outl)

p3 <- ggplot(data = meta_agg, mapping = aes(x = celltype, y = percentage, fill = clonestatus)) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.7) +
  facet_wrap(~ outcome,
             nrow = 1) +
  scale_fill_manual(values = bcrcol1) +
  labs(y = "Distribution of clone status (%)", x = "", title = "") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) + th
p3
