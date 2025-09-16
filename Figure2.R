library(Seurat)
library(rio)
library(tidyverse)
library(ggpubr)
library(ggplot2)
filep="./MiddleFile/"
sce.all=readRDS(paste0(filep,'4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))
outcomeCol=c('#7ac7e2','#f9d580','#e3716e')

th=theme(text = element_text(#family = "Arial",
  color="black"), 
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 0.5),  
  axis.ticks = element_line(linewidth =0.25) , 
  axis.text = element_text(size = 10, face = "bold", color = "black"))

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

outcomeCol=c('#7ac7e2','#f9d580','#e3716e')
## Figure 2A
# This approach was applied by Jaclyn M. L. Walsh in his manuscript "Variants and vaccines impact nasal immunity over three waves of SARS-CoV-2" (DOI: 10.1038/s41590-024-02052-z). Code from this paper (https://github.com/jo-m-lab/SARSCoV2-Variants-Vaccines-Paper) was adapted for the following analysis.
obj <-sce.all
obj$Detailed_Annotation=obj$celltype
obj$Participant=obj$orig.ident
obj$Coarse_Annotation=obj$celltype_major
obj$Outcome_Group=obj$outcome

obj@meta.data$tierNident <- obj@meta.data$celltype
obj@meta.data$CellID <- rownames(obj@meta.data)
obj@meta.data$sample <- obj@meta.data$orig.ident
meta=obj@meta.data

# First, make table of cells to calculate ratios. 
cells.md = obj@meta.data
cells.md$Participant = factor(as.character(cells.md$Participant), levels = sort(unique(as.character(cells.md$Participant))))
cells.md$Detailed_Annotation = factor(cells.md$Detailed_Annotation, levels = ctl)

#Calculate the frequencies and multiply by 500 (abundance)
labeled.freq.table = cells.md %>% select(Participant, Coarse_Annotation, Detailed_Annotation, Outcome_Group) %>% 
  group_by(Participant, Detailed_Annotation, .drop = FALSE) %>% dplyr::summarise(n = n()) %>%
  dplyr::mutate(abundance = n/sum(n)*500) %>%
  dplyr::mutate(log_abundance = log(abundance + 1)) 

#Create the cell type abundance matrices (22 annotated clusters, 72 Participants)
abundance.mtx = matrix(labeled.freq.table$abundance, byrow = FALSE, nrow = length(unique(obj$Detailed_Annotation)), 
                       dimnames = list(levels(labeled.freq.table$Detailed_Annotation), levels(labeled.freq.table$Participant)))

# Create associated metadata table
meta_sub <- meta[,c("Participant",'Outcome_Group')] 
meta_sub=meta_sub[!duplicated(meta_sub$Participant),]
rownames(meta_sub)=meta_sub$Participant

#Create combined abundance matrix
meta_sub <- meta_sub[match(colnames(abundance.mtx), meta_sub$Participant), ]
identical(meta_sub$Participant,colnames(abundance.mtx))
pt.abundance.df <- cbind(meta_sub[,c("Outcome_Group")], 
                         t(as.data.frame(abundance.mtx)))

#Run PCA on the log transformed matrix
abundance.mtx.log <- log(abundance.mtx + 1)
abundance.pca <- prcomp(t(abundance.mtx.log), center = TRUE, scale = TRUE)

var_exp <- abundance.pca$sdev^2 / sum(abundance.pca$sdev^2)
pc1_var <- round(var_exp[1] * 100, 1)
pc2_var <- round(var_exp[2] * 100, 1)

#Plot compositional PCA with ellipses 
meta_sub$Outcome_Group=factor(meta_sub$Outcome_Group,levels = c('Control','Alive','Deceased'))
p1=ggbiplot::ggbiplot(abundance.pca, var.axes = FALSE, groups = meta_sub$Outcome_Group, ellipse = TRUE, ellipse.fill = F) +  
  scale_color_manual(values = outcomeCol) + 
  xlab(paste0("Compositional PC1 (", pc1_var, "%)")) +
  ylab(paste0("Compositional PC2 (", pc2_var, "%)"))+
  theme_classic() + 
  theme(axis.title = element_text(face = "bold", color = "black", size = 12), 
        axis.text = element_text(color = "black", size = 10), 
        legend.title =element_blank() ,
        legend.text = element_text(color = "black", size = 10))+
  th
p1

## Figure 2C-D
ctl_big=c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono','Platelet')
ctc_big=c('#f091a0',"#d35b7e",'#fab378','#80b5b8','#caadd8','#80b1d2','#6B4C3D', '#b9181a')

pca.loading.df <- as.data.frame(abundance.pca$rotation ) %>%
  rownames_to_column("Detailed_Annotation") %>% 
  arrange(PC1)

#add cell type so bars can be colored by cell type 
ca=cells.md[!duplicated( cells.md$celltype),c('celltype','celltype_major')]
ca=ca[match(pca.loading.df$Detailed_Annotation,ca$celltype),]
ca$celltype=as.vector(ca$celltype)
identical(ca$celltype,pca.loading.df$Detailed_Annotation)
pca.loading.df$Coarse_Annotation <- ca$celltype_major

#set order of cell types 
pca.loading.df$Coarse_Annotation <- factor(pca.loading.df$Coarse_Annotation, levels = ctl_big)

#plot loadings for compositional pc1 
p2.1=ggplot(pca.loading.df, 
            aes(x = reorder(Detailed_Annotation, PC1), y = PC1, fill = Coarse_Annotation)) + 
  geom_bar(stat= "identity", width = 0.85) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_fill_manual(values = ctc_big) + 
  labs(x='',y='Cellular composition PC1 loading')+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, color = "black"), 
        axis.text.y = element_text(color ="black"), 
        axis.title.x = element_blank(),
        legend.position = 'top')+th

#plot loadings for compositional pc2 
p2.2=ggplot(pca.loading.df, 
            aes(x = reorder(Detailed_Annotation, PC2), y = PC2, fill = Coarse_Annotation)) + 
  geom_bar(stat= "identity", width = 0.85) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") + 
  scale_fill_manual(values =ctc_big) + 
  labs(x='',y='Cellular composition PC2 loading')+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, color = "black"), 
        axis.text.y = element_text(color = "black"), 
        axis.title.x = element_blank(),
        legend.position = 'top')+th

p2.1/p2.2

## Figure 2B Hierarchical clustering based on outcome group 
library(ggdendro)
abundance.pca.centroids <- aggregate(abundance.pca$x, list(Outcome_Group = meta_sub$Outcome_Group), mean) %>%
  tibble::column_to_rownames("Outcome_Group")

dist.conditions <- dist(abundance.pca.centroids)

hc <- hclust(dist.conditions)
hc.dendro <- dendro_data(hc)

p3=ggdendrogram(hc.dendro) + 
  geom_text(data = hc.dendro$labels, aes(x = x, y = y, label=label, color = label), vjust = 1.5) + 
  scale_color_manual(values =outcomeCol) + 
  scale_x_continuous(expand = c(.2,0)) +
  ylab("Height") + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.line.y = element_line(color = "gray20"),
        axis.ticks.y = element_line(color = "gray20"),
        legend.position = "none")
p3