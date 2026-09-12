library(Seurat)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(patchwork)

filep="./MiddleFile/"

sce.all=readRDS(paste0(filep,'R1_4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))
load(paste0(filep,"miloR_group_result.RData"))

th <- theme(text = element_text(family = "Arial"),
            axis.title = element_text(family = "Arial",size = 10,face = "bold"),
            plot.title = element_text(family = "Arial",size=12,face = "bold"),
            plot.caption = element_text(family = "Arial",size=8),
            axis.line = element_line(linewidth = 1),
            axis.ticks = element_line(linewidth = 0.5),
            axis.text = element_text(family = "Arial",size = 10, face = "bold", color = "black"))

ctl_big=c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono','Platelet')
ctc_big=c('#f091a0',"#d35b7e",'#fab378','#80b5b8','#caadd8','#80b1d2','#6B4C3D', '#b9181a')
sce.all$celltype_major=factor(sce.all$celltype_major,levels = ctl_big)

outcomeCol <- c('#7ac7e2','#f9d580','#e3716e')
names(outcomeCol) <- c("Control","Alive","Deceased")

## Figure 2A
plot_df <- sce.all@meta.data %>%
  dplyr::group_by(orig.ident, celltype_major, outcome) %>%
  dplyr::tally(name = "count") %>%
  dplyr::ungroup()

sample_order <- sce.all@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(outcome = first(outcome), .groups="drop") %>%
  dplyr::mutate(outcome = factor(outcome, levels = c("Control","Alive","Deceased"))) %>%
  dplyr::arrange(outcome)

plot_df$orig.ident <- factor(plot_df$orig.ident, levels = sample_order$orig.ident)

sample_color_map <- outcomeCol[as.character(sample_order$outcome)]
names(sample_color_map) <- sample_order$orig.ident

p2A <- ggplot(plot_df, aes(x = orig.ident, y = count, fill = celltype_major)) +
  geom_col(position = "fill", linewidth = 0.2) +
  scale_fill_manual(values = ctc_big) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Donor", y = "Proportion of Cells", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9, color = sample_color_map),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.box = "vertical"
  ) + th

print(p2A)

## Figure 2B
df=import(paste0(filep,'meta-total.xlsx'))
colnames(df)
table(df$celltype_major)
b=df[!duplicated(df$orig.ident), c('orig.ident','outcome')]
result <- df %>%
  group_by(orig.ident, celltype_major) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
data=merge(result,b,by="orig.ident")

plot_boxplot <- function(data, ct) {
  cell_data <- data[data$celltype_major == ct, ]
  max_y <- max(cell_data$Proportion, na.rm = TRUE)
  
  p <- ggplot(cell_data, aes(x = outcome, y = Proportion, fill = outcome)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2, size = 1, alpha = 0.7) + 
    theme_classic() +  
    theme(
      axis.title.x = element_blank(),  
      plot.title = element_text(hjust = 0.5) 
    ) + 
    labs(title = ct, y = "Proportion", x = '') +  
    scale_fill_manual(values = c("Control" = '#7ac7e2', "Alive" = '#f9d580', "Deceased" = '#e3716e')) +  
    geom_signif(
      comparisons = list(
        c("Control", "Alive"),
        c("Alive", "Deceased"),
        c("Control", "Deceased")
      ),
      map_signif_level = TRUE,
      y_position = c(max_y * 1.1, max_y * 1.1, max_y * 1.2)  
    ) + 
    th
}

data$outcome=factor(data$outcome,levels = c('Control','Alive','Deceased'))
ct <- c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono','Platelet')
plots <- lapply(ct, function(cyt) plot_boxplot(data, cyt))
p2B <- wrap_plots(plots, ncol = 4)
p2B

## Figure 2D
color_control  <- "#7ac7e2"
color_alive    <- "#fabb6e"
color_deceased <- "#b7282e"
color_grey     <- "gray80"

plot_beeswarm <- function(da_result, title, type) {

  da_plot <- da_result

  da_plot$celltype <- factor(
    da_plot$celltype,
    levels = rev(c('B_c01-mature naive B','B_c02-memory B',
                   'B_c03-active Plasma','B_c04-terminal Plasma',
                   'T_CD4_c01-LEF1','T_CD4_c02-AQP3','T_CD4_c03-FOS','T_CD4_c04-FOXP3',
                   'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK','T_CD8_c03-TYROBP','T_CD8_c04-ZNF683','T_CD8_c05-RORA' ,
                   'NK_c01-CD56high','NK_c02-CD56low',
                   'DC-mDC',
                   'Mono_c01-Classical','Mono_c02-Intermediate','Mono_c03-Non_Classical',
                   'Platelet'
    )
  ))

  if(type == "Alive_vs_Control"){
    da_plot <- da_plot %>% mutate(
      point_color = case_when(
        FDR < 0.05 & logFC < -1 ~ "Enriched_in_Control",
        FDR < 0.05 & logFC > 1 ~ "Enriched_in_Alive",
        TRUE ~ "Not_significant"
      ),
      bar_color = case_when(
        (any(FDR < 0.05 & logFC < -1)) & mean(logFC) < 0 ~ "Enriched_in_Control",
        (any(FDR < 0.05 & logFC > 1)) & mean(logFC) > 0 ~ "Enriched_in_Alive",
        TRUE ~ "Neutral"
      ), .by = celltype
    )
    col_values <- c(
      Enriched_in_Control = color_control,
      Enriched_in_Alive = color_alive,
      Not_significant = color_grey
    )
    bar_values <- c(
      Enriched_in_Control = color_control,
      Enriched_in_Alive = color_alive,
      Neutral = "white"
    )
  }

  else if(type == "Deceased_vs_Control"){
    da_plot <- da_plot %>% mutate(
      point_color = case_when(
        FDR < 0.05 & logFC < -1 ~ "Enriched_in_Control",
        FDR < 0.05 & logFC > 1 ~ "Enriched_in_Deceased",
        TRUE ~ "Not_significant"
      ),
      bar_color = case_when(
        (any(FDR < 0.05 & logFC < -1)) & mean(logFC) < 0 ~ "Enriched_in_Control",
        (any(FDR < 0.05 & logFC > 1)) & mean(logFC) > 0 ~ "Enriched_in_Deceased",
        TRUE ~ "Neutral"
      ), .by = celltype
    )
    col_values <- c(
      Enriched_in_Control = color_control,
      Enriched_in_Deceased = color_deceased,
      Not_significant = color_grey
    )
    bar_values <- c(
      Enriched_in_Control = color_control,
      Enriched_in_Deceased = color_deceased,
      Neutral = "white"
    )
  }

  else if(type == "Deceased_vs_Alive"){
    da_plot <- da_plot %>% mutate(
      point_color = case_when(
        FDR < 0.05 & logFC < -1 ~ "Enriched_in_Alive",
        FDR < 0.05 & logFC > 1 ~ "Enriched_in_Deceased",
        TRUE ~ "Not_significant"
      ),
      bar_color = case_when(
        (any(FDR < 0.05 & logFC < -1)) & mean(logFC) < 0 ~ "Enriched_in_Alive",
        (any(FDR < 0.05 & logFC > 1)) & mean(logFC) > 0 ~ "Enriched_in_Deceased",
        TRUE ~ "Neutral"
      ), .by = celltype
    )
    col_values <- c(
      Enriched_in_Alive = color_alive,
      Enriched_in_Deceased = color_deceased,
      Not_significant = color_grey
    )
    bar_values <- c(
      Enriched_in_Alive = color_alive,
      Enriched_in_Deceased = color_deceased,
      Neutral = "white"
    )
  }

  p <- ggplot(da_plot, aes(x = logFC, y = celltype)) +
    geom_quasirandom(
      data = subset(da_plot, point_color == "Not_significant"),
      color = color_grey, size = 0.9, alpha = 0.6
    ) +
    geom_quasirandom(
      data = subset(da_plot, point_color != "Not_significant"),
      aes(color = point_color), size = 1.1, alpha = 0.8
    ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-2, 2), color = "black", linetype =  "dashed", linewidth = 0.6, alpha=0.7) +

    geom_tile(aes(x = -7, fill = bar_color), width = 0.6, height = 0.4) +

    scale_color_manual(values = col_values, guide = "none") +
    scale_fill_manual(values = bar_values, guide = "none") +
    scale_x_continuous(limits = c(-7.5, 5.5), breaks = seq(-5, 5, by=2.5)) +
    labs(x = "Log fold change", y = "", title = title) +
    theme_bw() +
    theme(axis.title = element_text(size = 10,face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          plot.caption = element_text(size=8),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 0.5),
          axis.text = element_text(size = 10, face = "bold", color = "black"),
          panel.grid = element_blank(),
          plot.margin = margin(10,10,10,40)
          )

  return(p)
}

p1 <- plot_beeswarm(da_Alive_vs_Control,    "Control vs Alive\nDA direction",        "Alive_vs_Control")
p2 <- plot_beeswarm(da_Deceased_vs_Control, "Control vs Deceased\nDA direction",     "Deceased_vs_Control")
p3 <- plot_beeswarm(da_Deceased_vs_Alive,   "Alive vs Deceased\nDA direction",       "Deceased_vs_Alive")

p_combined <- p1 + p2 + p3 + plot_layout(nrow = 1)

print(p_combined)
