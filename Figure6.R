library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(rio)
library(pheatmap)
library(ggrepel)
library(scales)

th=theme(text = element_text(
  color="black"),
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1.5),
  axis.ticks = element_line(linewidth = 1),
  axis.text = element_text(size = 10, face = "bold", color = "black"))

colt=c(
  '#fcedbe','#fab378','#f49512','#ffd377',
  '#6da5ad','#d9e773',
  '#ccd4b5','#889851',
  '#BAD8D2'
)
outcomeCol=c('#7ac7e2','#f9d580','#e3716e')
outl=c('Control','Alive','Deceased')
names(outcomeCol)=outl
hmc=c("#aceefe", "white",'#ffd377')
filep="./MiddleFile/"

ct_tl=c('T_CD4_c01-LEF1','T_CD4_c02-AQP3','T_CD4_c03-FOS','T_CD4_c04-FOXP3',
        'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK','T_CD8_c03-TYROBP','T_CD8_c04-ZNF683','T_CD8_c05-RORA')

## Figure 6A
data=import(paste0(filep,'Fig6A_Tcell_function_ternary_SourceData.xlsx'),which='Plot_data')
data$T_celltype=factor(data$T_celltype,levels=ct_tl)

angles <- (c(0, 120, 240) * pi / 180)

axis_df <- data.frame(
  x_start = rep(0, 3),
  y_start = rep(0, 3),
  x_end = c(1, cos(angles[2]), cos(angles[3])),
  y_end = c(0, sin(angles[2]), sin(angles[3]))
)

axis_labels <- data.frame(
  x = c(1.1, 1.1*cos(angles[2]), 1.1*cos(angles[3])),
  y = c(0, 1.1*sin(angles[2]), 1.1*sin(angles[3])),
  label = c("Cytotoxic", "Naive", "Regulatory")
)

ggplot() +
  geom_segment(
    data = axis_df,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth = 0.8,
    color = "gray30"
  ) +
  geom_text(
    data = axis_labels,
    aes(x, y, label = label),
    size = 5,
    color = "black"
  ) +
  geom_point(
    data = data,
    aes(x, y, color = T_celltype),
    size = 3
  ) +
  geom_text_repel(
    data = data,
    aes(x, y, label = T_celltype, color = T_celltype),
    size = 4,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  th+
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = 'none'
  ) +
  scale_color_manual(values = colt)

## Figure 6B
sce.all=readRDS(paste0(filep,'R1_4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))
subtypes <- ct_tl
meta <- sce.all@meta.data

sample_total <- meta %>%
  group_by(orig.ident, outcome) %>%
  summarise(total_cells = n(), .groups = "drop")

prop <- meta %>%
  filter(celltype %in% subtypes) %>%
  group_by(orig.ident, outcome, celltype) %>%
  summarise(n_t = n(), .groups = "drop") %>%
  right_join(
    expand_grid(
      orig.ident = unique(meta$orig.ident),
      celltype = subtypes
    ),
    by = c("orig.ident", "celltype")
  ) %>%
  left_join(sample_total, by = "orig.ident") %>%
  mutate(
    outcome = coalesce(outcome.x, outcome.y),
    n_t = replace_na(n_t, 0),
    prop = n_t / total_cells
  ) %>%
  select(orig.ident, outcome, celltype, n_t, total_cells, prop)

prop$outcome <- factor(prop$outcome, levels = outl)

ggplot(prop, aes(x = outcome, y = prop, fill = outcome)) +
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
  facet_wrap(~celltype, scales = "free_y",ncol = 5) +
  scale_fill_manual(values = outcomeCol) +
  scale_color_manual(values = outcomeCol) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    x = "",
    y = "Proportion of all cells",
    title = ""
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

## Figure 6C
naive_subsets <- c(
  "T_CD4_c01-LEF1",
  "T_CD4_c02-AQP3",
  "T_CD8_c01-LEF1"
)

donor_score <- read.csv(paste0(filep, "NaiveT_functional_scores_donor_level.csv"))
dunn_res    <- read.csv(paste0(filep, "NaiveT_scores_stats_Dunn.csv"))

scores_show <- c("Apoptosis", "Quiescence")

plot_df <- donor_score %>%
  filter(celltype %in% naive_subsets) %>%
  select(orig.ident, outcome, celltype, all_of(scores_show)) %>%
  pivot_longer(
    cols = all_of(scores_show),
    names_to = "score",
    values_to = "value"
  ) %>%
  mutate(
    celltype = factor(celltype, levels = naive_subsets),
    outcome  = factor(outcome, levels = outl),
    score    = factor(score, levels = scores_show)
  )

make_dumbbell <- function(score_name, df_all, dunn_all) {

  sub_df <- df_all %>% filter(score == score_name)

  med_df <- sub_df %>%
    group_by(celltype, outcome) %>%
    summarise(med = median(value), .groups = "drop") %>%
    mutate(x = as.numeric(outcome))

  sig <- dunn_all %>%
    filter(score == score_name, p.adj < 0.05) %>%
    mutate(
      xmin = pmin(match(group1, outl), match(group2, outl)),
      xmax = pmax(match(group1, outl), match(group2, outl)),
      label = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        TRUE          ~ "*"
      ),
      celltype = factor(celltype, levels = naive_subsets)
    )

  if (nrow(sig) > 0) {
    yr <- sub_df %>%
      group_by(celltype) %>%
      summarise(ymax = max(value), ymin = min(value), .groups = "drop") %>%
      mutate(yoff = 0.08 * (ymax - ymin))

    sig <- sig %>%
      left_join(yr, by = "celltype") %>%
      group_by(celltype) %>%
      arrange(desc(xmax - xmin), .by_group = TRUE) %>%
      mutate(ypos = ymax + yoff * row_number()) %>%
      ungroup()
  }

  p <- ggplot(sub_df, aes(x = outcome, y = value)) +
    geom_jitter(aes(color = outcome), width = 0.18, size = 0.9, alpha = 0.55) +
    geom_line(
      data = med_df,
      aes(x = x, y = med, group = celltype),
      color = "grey50", linewidth = 0.6
    ) +
    geom_point(
      data = med_df,
      aes(x = x, y = med, color = outcome),
      size = 3
    ) +
    scale_color_manual(values = outcomeCol) +
    facet_wrap(~ celltype, nrow = 1) +
    theme_classic(base_size = 11) +
    labs(x = NULL, y = NULL, title = score_name) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey92", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black")
    )

  if (nrow(sig) > 0) {
    p <- p +
      geom_segment(
        data = sig,
        aes(x = xmin, xend = xmax, y = ypos, yend = ypos),
        inherit.aes = FALSE,
        linewidth = 0.5, color = "black"
      ) +
      geom_segment(
        data = sig,
        aes(x = xmin, xend = xmin, y = ypos, yend = ypos - 0.3 * (yoff)),
        inherit.aes = FALSE,
        linewidth = 0.5, color = "black"
      ) +
      geom_segment(
        data = sig,
        aes(x = xmax, xend = xmax, y = ypos, yend = ypos - 0.3 * (yoff)),
        inherit.aes = FALSE,
        linewidth = 0.5, color = "black"
      ) +
      geom_text(
        data = sig,
        aes(x = (xmin + xmax) / 2, y = ypos + 0.4 * yoff, label = label),
        inherit.aes = FALSE,
        size = 4, color = "black"
      ) +
      expand_limits(y = max(sig$ypos) * 1.02)
  }

  return(p)
}

p_apop <- make_dumbbell("Apoptosis",  plot_df, dunn_res)

## Figure 6D
p_quie <- make_dumbbell("Quiescence", plot_df, dunn_res)

p_dumb <- (p_apop | p_quie)
p_dumb

## Figure 6E
meta=import(paste0(filep,'CD4_8-T-score_meta.xlsx'))
meta1=meta[meta$celltype %in% c('T_CD4_c03-FOS'),]

selected_columns <-c("TCR signaling","Costimulatory molecules",
                     "Activation/Effector function",
                     'Cytotoxicity',
                     'Cytokine/Cytokine receptor'
)

mean_values <- meta1 %>%
  group_by(outcome) %>%
  summarise(across(all_of(selected_columns), mean, na.rm = TRUE))
mean_values=as.data.frame(mean_values)

rownames(mean_values)=mean_values$outcome
mean_values=mean_values[,-1]
hm=t(mean_values)
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = 'T_CD4_c03-FOS')

## Figure 6F
meta=import(paste0(filep,'CD8_T-score_meta.xlsx'))
selected_columns <-c(
  'TCR Signaling',
  "NFKB Signaling",
  'Cytotoxicity',"Stress response",
  'Exhaustion'
)

cn=c('T_CD8_c02-GZMK')
meta1=meta[meta$celltype %in% cn,]
mean_values <- meta1 %>%
  group_by(outcome) %>%
  summarise(across(all_of(selected_columns), mean, na.rm = TRUE))
mean_values=as.data.frame(mean_values)
rownames(mean_values)=mean_values$outcome
mean_values=mean_values[,-1]
hm=t(mean_values)
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = cn)

cn=c('T_CD8_c03-TYROBP')
meta1=meta[meta$celltype %in% cn,]
mean_values <- meta1 %>%
  group_by(outcome) %>%
  summarise(across(all_of(selected_columns), mean, na.rm = TRUE))
mean_values=as.data.frame(mean_values)
rownames(mean_values)=mean_values$outcome
mean_values=mean_values[,-1]
hm=t(mean_values)
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = cn)

cn=c('T_CD8_c04-ZNF683')
meta1=meta[meta$celltype %in% cn,]
mean_values <- meta1 %>%
  group_by(outcome) %>%
  summarise(across(all_of(selected_columns), mean, na.rm = TRUE))
mean_values=as.data.frame(mean_values)
rownames(mean_values)=mean_values$outcome
mean_values=mean_values[,-1]
hm=t(mean_values)
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = cn)

cn=c('T_CD8_c05-RORA')
meta1=meta[meta$celltype %in% cn,]
mean_values <- meta1 %>%
  group_by(outcome) %>%
  summarise(across(all_of(selected_columns), mean, na.rm = TRUE))
mean_values=as.data.frame(mean_values)
rownames(mean_values)=mean_values$outcome
mean_values=mean_values[,-1]
hm=t(mean_values)
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = cn)

## Figure 6H
meta=import(paste0(filep,'CD4_c04-ImmuneSuppressioin-score_meta.xlsx'))
meta1=meta[meta$celltype %in% c('T_CD4_c04-FOXP3'),]
meta2=meta[!duplicated(meta$orig.ident),]
meta2=meta2[,c('orig.ident','outcome')]

selected_columns <-c('Immune_suppression')

mean_values <- meta1 %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(selected_columns), mean, na.rm = TRUE))
mean_values=as.data.frame(mean_values)

df=merge(mean_values,meta2,by='orig.ident')
df$outcome=factor(df$outcome,levels = c('Control','Alive','Deceased'))

ggplot(df, aes(x = outcome, y =Immune_suppression, fill =outcome)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Immune suppression",title = 'T_CD4_c04-FOXP3') +
  theme_classic() +
  theme(legend.position = "none")+th

## Figure 6I
df=rio::import(paste0(filep,'TCR_clonalDiversity-sample.xlsx'))
meta=sce.all@meta.data
meta1=meta[!duplicated(meta$orig.ident) ,c('orig.ident','outcome')]
df1=merge(df[,c('Shannon','Inv.Simpson','Sample')],meta1,by.x='Sample',by.y='orig.ident')

pt1.1=ggboxplot(data = df1,
                x = 'outcome',
                y = 'Shannon',
                color = 'gray40',
                fill = NA,
                outlier.shape = NA
) +
  geom_jitter(aes(x = outcome, y = Shannon, color = outcome),
              width = 0.2, size = 2, alpha = 1) +
  scale_color_manual(values = c('Control' = '#7ac7e2',
                                'Alive' = '#f9d580',
                                'Deceased' = '#e3716e')) +
  labs(title = '', x = '', y = 'Shannon index') +
  theme_classic()+
  th

pt1.2=ggboxplot(data = df1,
                x = 'outcome',
                y = 'Inv.Simpson' ,
                color = 'gray40',
                fill = NA,
                outlier.shape = NA
) +
  geom_jitter(aes(x = outcome, y = Inv.Simpson, color = outcome),
              width = 0.2, size = 2, alpha = 1) +
  scale_color_manual(values = c('Control' = '#7ac7e2',
                                'Alive' = '#f9d580',
                                'Deceased' = '#e3716e')) +
  labs(title = '', x = '', y = 'Inv.Simpson index') +
  theme_classic()+
  th


## Figure 6J
meta=import(paste0(filep,'meta-T_TCRmerged.xlsx'))
tcrcol=c("#cbe3e6",
         '#889851','#d9e773','#fcedbe','#ffd377','#fab378')

cell_types <- meta
cell_types$outcome=factor(cell_types$outcome,levels=c("Control","Alive","Deceased"))
cell_types=cell_types[cell_types$cloneSize !="None",]
cell_types$size=ifelse(cell_types$clonalFrequency==0,"Non-VDJ",ifelse(cell_types$clonalFrequency==1,"unique",
                       ifelse(cell_types$clonalFrequency <= 5, "2-5",ifelse(cell_types$clonalFrequency<=20,"6-20",
                                                                            ifelse(cell_types$clonalFrequency<=50,"21-50",
                                                                                   ifelse(cell_types$clonalFrequency<=100,"51-100",
                                                                                          ">100")))

                       )))

cell_types$size=factor(cell_types$size,
                       levels=c("unique","2-5",'6-20',
                                '21-50','51-100','>100'))
cell_types_agg <- cell_types %>%
  group_by(outcome, size) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(outcome) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()
cell_types_agg$outcome=factor(cell_types_agg$outcome,levels = c("Control","Alive","Deceased"))

ggplot(data = cell_types_agg, mapping = aes(x = "", y = prop, fill = size)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = tcrcol) +
  labs(y = "", x = "", title = "") +
  theme_classic() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 9,face = 'bold'),
        strip.background = element_blank()) +
  geom_label_repel(aes(label = scales::percent(prop, accuracy = 0.01)),
                   position = position_stack(vjust = 0.5),
                   box.padding = 0.5,
                   size = 3,
                   show.legend = FALSE) +
  facet_wrap(~ outcome, nrow = 1)
