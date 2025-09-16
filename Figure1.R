# Figure 1B heatmap
library(ComplexHeatmap)
library(circlize)
library(rio)
library(scales)

a=import_list("PBMC-meta information.xlsx")[[1]]
colnames(a)
data <- data.frame(
  Sample=a$SampleID,
  Group = a$Group,
  Outcome = a$Outcome,
  Sex = a$Sex,
  Age = a$Age_group
)

data_t <- t(data[,-1])
colnames(data_t) <- data$Sample

group_col <- c("Healthy Control" = "#C2D7F3", "Severe Pneumonia" = "#F7CF83")
outcome_col <- c("Control" ="#C2D7F3", "Alive" = '#fabb6e', "Deceased" = '#cc625f')
sex_col <- c("Female" = "#f8cbe0", "Male" = "#97ceff")
show_col(sex_col)
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
    Group = factor(data$Group, levels = c("Healthy Control", "Severe Pneumonia")),
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

Heatmap(as.matrix(data_t), name = "Heatmap", 
          top_annotation = ha, 
          column_names_side = "top",
          show_row_names = TRUE, show_column_names = TRUE,
          column_names_gp = gpar(fontsize = 10), 
          row_names_gp = gpar(fontsize = 10)
)  


# Figure 1C
library(Seurat)
library(ggplot2)
th=theme(text = element_text(#family = "Arial",
  color="black"),  
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1),  
  axis.ticks = element_line(linewidth = 0.5) , 
  axis.text = element_text(size = 10, face = "bold", color = "black"))

filep="./MiddleFile/"

sce.all=readRDS(paste0(filep,'4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))

table(sce.all$celltype)
ctl=c('B_c01-mature naive B','B_c02-memory B',
      'B_c03-Plasmablast',
      'B_c04-active Plasma','B_c05-terminal Plasma',
      'T_CD4_c01-LEF1','T_CD4_c02-AQP3','T_CD4_c03-FOS','T_CD4_c04-FOXP3',
      'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK','T_CD8_c03-TYROBP','T_CD8_c04-ZNF683','T_CD8_c05-RORA' ,
      'NK_c01-CD56high','NK_c02-CD56low',
      'DC-mDC',
      'Mono_c01-CD14-VCAN','Mono_c02-CD14-CCL3','Mono_c03-CD14-HLA-DPB1','Mono_c04-CD14-CD16',
      'Platelet'
)
sce.all$celltype=as.vector(sce.all$celltype)
sce.all$celltype=factor(sce.all$celltype, levels = ctl)

ctc=c( # B & Plasma
  '#f9d3e3','#cb7395','#ffadbb','#8f476d',"#d35b7e",
           #CD4+T
           '#fcedbe','#fab378','#f49512','#ffd377',
           #CD8+T
           '#6da5ad','#d9e773', '#ccd4b5','#889851', '#BAD8D2',
           #NK
           "#c9d4f7",'#caadd8',
           #DC
           '#b1ebfe',
           #Mono        
           '#d7cdbe','#6B4C3D','#cd9e80',
           '#997E68',
           #Platelet
           '#b9181a')
names(ctc)=ctl

DimPlot(sce.all,group.by = "celltype",reduction = "umap_harmony",
        cols = ctc) +
  labs(x='UMAP1',y='UMAP2',title = '')+th+
  theme(legend.position = 'bottom', legend.justification = 'center') +
  guides(color = guide_legend(ncol = 2))  


## Figure 1D
library(dplyr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggbiplot)
library(ggdendro)
library(dunn.test)
library(effsize)
library(ARBOL)

ctl_big=c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono','Platelet')
ctc_big=c('#f091a0',"#d35b7e",'#fab378','#80b5b8','#caadd8','#80b1d2','#6B4C3D', '#b9181a')
sce.all$celltype_major=factor(sce.all$celltype_major,levels = ctl_big)

pal.celltypes <-ctc_big
names(pal.celltypes) <- ctl_big
pal.celltypes

obj <-sce.all

obj$Detailed_Annotation=obj$celltype
obj$Participant=obj$orig.ident
obj$Coarse_Annotation=obj$celltype_major
obj$Outcome_Group=obj$outcome

obj@meta.data$tierNident <- obj@meta.data$celltype
obj@meta.data$CellID <- rownames(obj@meta.data)
obj@meta.data$sample <- obj@meta.data$orig.ident
meta <- obj@meta.data

obj <- ARBOLcentroidTaxonomy(obj,
                             categories = c('Participant', 'Outcome_Group', 'Coarse_Annotation'),
                             diversities = c('Participant', 'Coarse_Annotation'),
                             counts = c('Participant', 'Outcome_Group', 'Coarse_Annotation'),
                             tree_reduction='harmony', 
                             centroid_method = 'median', 
                             distance_method ="euclidean",
                             hclust_method='ward.D2',  #'complete',
                             nboot=1)

ggraph(obj@misc$tax_ggraph, layout='dendrogram') + 
  geom_edge_elbow2(aes(color = node.Coarse_Annotation_majority), width = 2) +
  scale_edge_colour_manual(values = pal.celltypes) + 
  guides(edge_colour = "none") + 
  theme_void() + 
  expand_limits(y = -5) + 
  new_scale_color() +
  guides(color = FALSE) +  
  geom_node_text(aes(filter = leaf, label = name, color = Coarse_Annotation_majority),
                 nudge_y = -1.2, vjust = 0.5, hjust = 0, angle = 270, size = 5) +
  scale_color_manual(values = pal.celltypes) +
  new_scale_color() +
  geom_node_text(aes(filter = leaf, label = n),
                 color = 'grey30', nudge_y = -0.3, vjust = 0.5, hjust = 0, size = 5, angle = 270) +
  new_scale_color() +
  geom_node_point(aes(filter = leaf, color = Participant_diversity),
                  size = 5, shape = 'square') +
  scale_color_gradient(
    low = 'grey90', high = 'grey10',
    limits = c(0, 1), breaks = seq(0, 1, 0.1) 
  ) +
  
  theme(
    legend.text = element_text(size = 15), 
    legend.title = element_text(size = 15, face = "bold")
  ) + 
  expand_limits(y = -10)
