library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(rio)
library(ggraph)
library(ggnewscale)

filep="./MiddleFile/"

sce.all=readRDS(paste0(filep,'R1_4.publicHC_PBMC-sce.all_qc_cluster-join-celltype-2.rds'))
a=import_list(paste0(filep,'PBMC-meta information-250423.xlsx'))[[1]]
arbol=readRDS(paste0(filep,'R1_Fig1_ARBOL_tree_result.rds'))

th=theme(text = element_text(color="black"),
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1),
  axis.ticks = element_line(linewidth = 0.5),
  axis.text = element_text(size = 10, face = "bold", color = "black"))

ctl=c('B_c01-mature naive B','B_c02-memory B',
      'B_c03-active Plasma','B_c04-terminal Plasma',
      'T_CD4_c01-LEF1','T_CD4_c02-AQP3','T_CD4_c03-FOS','T_CD4_c04-FOXP3',
      'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK','T_CD8_c03-TYROBP','T_CD8_c04-ZNF683','T_CD8_c05-RORA' ,
      'NK_c01-CD56high','NK_c02-CD56low',
      'DC-mDC',
      'Mono_c01-Classical','Mono_c02-Intermediate','Mono_c03-Non_Classical',
      'Platelet'
)
sce.all$celltype=as.vector(sce.all$celltype)
sce.all$celltype=factor(sce.all$celltype, levels = ctl)

ctc=c('#f9d3e3','#cb7395'
           ,'#8f476d',"#d35b7e",
           '#fcedbe','#fab378','#f49512','#ffd377',
           '#6da5ad','#d9e773',
           '#ccd4b5','#889851',
           '#BAD8D2',
           "#c9d4f7",'#caadd8',
           '#b1ebfe',
           '#d7cdbe','#6B4C3D','#cd9e80',
           '#b9181a')
names(ctc)=ctl

ctl_big=c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono','Platelet')
ctc_big=c('#f091a0',"#d35b7e",'#fab378','#80b5b8','#caadd8','#80b1d2','#6B4C3D', '#b9181a')
sce.all$celltype_major=factor(sce.all$celltype_major,levels = ctl_big)

## Figure 1B
data <- data.frame(
  Sample=a$SampleID,
  Group = a$Sampletype,
  Outcome = a$Outcome,
  Sex = a$Sex,
  Age = a$Age_group
)

data_t <- t(data[,-1])
colnames(data_t) <- data$Sample

group_col <- c("scRNA-seq" ='#e68b81', "Flow cytometry" = '#f7df87')
outcome_col <- c("Control" ="#C2D7F3", "Alive" = '#fabb6e', "Deceased" = '#cc625f')
sex_col <- c("Female" = "#f8cbe0", "Male" = "#97ceff")
age_col <- c("20-29" = "#FAEE85",
             "30-39" = "#FED477",
             "40-49" = "#E4C455",
             "50-59" = "#eca680",
             "60-69" = "#F8C9D5",
             "70-79" = "#D8B8D6",
             "80-89" = "#EE84A8",
             "90-95" ='#8f476d' )

ha <- HeatmapAnnotation(
  df = data.frame(
    Group = factor(data$Group, levels = c("scRNA-seq", "Flow cytometry" )),
    Outcome = factor(data$Outcome, levels = c("Control", "Alive", "Deceased")),
    Sex = factor(data$Sex, levels = c("Female", "Male")),
    Age = factor(data$Age, levels = c("20-29","30-39","40-49","50-59", "60-69", "70-79", "80-89", "90-95"))
  ),
  col = list(
    Group = group_col,
    Outcome = outcome_col,
    Sex = sex_col,
    Age = age_col
  )
)

p=Heatmap(as.matrix(data_t), name = "Heatmap",
          top_annotation = ha,
          column_names_side = "top",
          show_row_names = TRUE, show_column_names = TRUE,
          column_names_gp = gpar(fontsize = 10),
          row_names_gp = gpar(fontsize = 10)
)
print(p)

## Figure 1C
DimPlot(sce.all,group.by = "celltype",reduction = "umap_harmony",
        raster=TRUE,
        cols = ctc) +
  labs(x='UMAP1',y='UMAP2',title = '')+th+
  theme(legend.position = 'bottom', legend.justification = 'center') +
  guides(color = guide_legend(ncol = 2))

## Figure 1D
pal.celltypes <-ctc_big
names(pal.celltypes) <- ctl_big

ggraph(arbol$tax_ggraph, layout='dendrogram') +
  geom_edge_elbow2(
    aes(color = as.character(node.Coarse_Annotation_majority)),
    width = 2
  ) +
  scale_edge_colour_manual(values = pal.celltypes, na.value = "grey50") +
  guides(edge_colour = "none") +

  theme_void() +
  expand_limits(y = -5) +
  new_scale_color() +

  geom_node_text(
    aes(filter = leaf, label = name, color = as.character(Coarse_Annotation_majority)),
    nudge_y = -1.2, vjust = 0.5, hjust = 0, angle = 270, size = 5
  ) +
  scale_color_manual(values = pal.celltypes, na.value = "grey50") +
  guides(color = "none") +
  new_scale_color() +

  geom_node_text(
    aes(filter = leaf, label = n),
    color = 'grey30', nudge_y = -0.3, vjust = 0.5, hjust = 0, angle = 270, size = 5
  ) +
  new_scale_color() +

  geom_node_point(
    aes(filter = leaf, fill = Participant_diversity),
    size = 5, shape = 22
  ) +
  scale_fill_gradient(
    low = 'grey90', high = 'grey10',
    limits = c(0.7, 1), breaks = seq(0.7, 1, 0.1)
  ) +

  theme(
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15, face = "bold")
  ) +
  expand_limits(y = -10)
