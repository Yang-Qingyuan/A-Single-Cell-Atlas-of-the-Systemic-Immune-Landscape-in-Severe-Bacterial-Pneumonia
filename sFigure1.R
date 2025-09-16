## sFigure 1A
library(rio)
library(ggplot2)
library(tidyverse)
th=theme(text = element_text(family = "Arial",color="black"),  
         axis.title = element_text(size = 10,face = "bold"),
         plot.title = element_text(size=12,face = "bold"),
         plot.caption = element_text(size=8),
         axis.line = element_line(size = 1), 
         axis.ticks = element_line(linewidth = 0.5) , 
         axis.text = element_text(size = 10, face = "bold", color = "black")
)

a=import_list("PBMC-meta information-250423.xlsx")[[1]]
outcomeCol=c('#7ac7e2','#f9d580','#e3716e')

df <- a %>% arrange(Age)
df$SampleID=factor(df$SampleID,levels = df$SampleID)
sp1.1= ggplot(df, aes(x = SampleID, y = Age, color = Outcome)) +
  geom_segment(aes(x = SampleID, xend = SampleID, y = 0, yend = Age, color = Outcome), size = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = outcomeCol, breaks = c("Control", "Alive", "Deceased")) +
  labs(x = "Patients ranked by age", y = "Age") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) + 
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank() 
  ) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "grey") + 
  guides(color = guide_legend(title = NULL))+ 
  th
sp1.1  

## sFigure 1B
library(tidyverse)
outl=c('Control','Alive','Deceased')
sexcol=c( "#f8cbe0",  "#97ceff")

count_data <- a %>% 
  group_by(Outcome, Sex) %>% 
  summarise(count = n()) %>% 
  ungroup()

count_data$Outcome=factor(count_data$Outcome, levels=outl)

chisq_test <- chisq.test(table(a$Outcome, a$Sex))
chisq_test$expected
p_value <- chisq_test$p.value

sp1.2=ggplot(count_data, aes(x =Outcome, y = count, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = sexcol) +
  labs(x = "", y = "Number  of Donors", fill = "Sex") +
  theme_classic()+
  th+
  theme(legend.position = "top")+
  guides(fill = guide_legend(title = NULL)) + 
  annotate("text", x = 0.5, y = max(count_data$count) + 1, label = paste("P-value:", round(p_value, 3)), size = 3, hjust = 0)
sp1.2

## sFigure 1C
library(Seurat)
filep="./MiddleFile/"
sce.all=readRDS(paste0(filep,'4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))

meta=sce.all@meta.data
meta1=meta
colnames(meta1)
meta1$outcome=factor(meta1$outcome,levels = outl)

sp1.3=ggplot(meta1, aes(x =outcome, y = nCount_RNA, fill = outcome)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Number of Counts") +
  theme_classic() +
  theme(legend.position = "none")+th
sp1.3

## sFigure 1D
colnames(meta1)
sp1.4=ggplot(meta1, aes(x = outcome, y = nFeature_RNA, fill =outcome)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Number of Genes") +
  theme_classic() +
  theme(legend.position = "none")+th
sp1.4

## sFigure 1E-H
Dot1c=c("#F8D941","#C73866")
Idents(sce.all)=sce.all$celltype_major

# B
sce=subset(sce.all,idents = c("B","Plasma"))
Idents(sce)=sce$celltype
sce$celltype=factor(sce$celltype,levels=c('B_c01-mature naive B','B_c02-memory B',
                                          'B_c03-Plasmablast',
                                          'B_c04-active Plasma','B_c05-terminal Plasma'))
genes = c(
  'CD19', 'CD79A','CD79B', 'MS4A1',# B
  # mature naive-B
  'IGHD','IGHM','CCR7','TCL1A',
  # Memory B
  'EBI3','FCRL4','DUSP4','GPR183',
  # Plasmablast
  "CXCR3", "IRF4", "KLF4",  'MKI67','PCNA',
  #Plasma active
  'PRDM1','XBP1','CCR10','TXNDC5',
  #Plasma terminal
  'LGALS3','CCL2','CCL5','ANXA1'
)

sp1.b <- DotPlot(sce, 
                 features =unique(genes),
                 cols =Dot1c,
                 assay='RNA'  )  + coord_flip()+
  labs(x="",y="",title = "")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "B cells")+
  th
sp1.b

# T
Idents(sce.all)=sce.all$celltype_major
sce=subset(sce.all,idents=c('CD4+T','CD8+T'))
sce$celltype=as.vector(sce$celltype)
sce$celltype <- factor(sce$celltype, levels = c(
  'T_CD4_c01-LEF1','T_CD4_c02-AQP3','T_CD4_c03-FOS','T_CD4_c04-FOXP3',
  'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK','T_CD8_c03-TYROBP','T_CD8_c04-ZNF683','T_CD8_c05-RORA' 
))
genes= c(
  'CD3D', 'CD3E',
  'CD4',
  'LEF1','AQP3','FOS','FOXP3',
  'CD8A','CD8B',
  'GZMK','TYROBP','ZNF683','RORA'
)
Idents(sce)=sce$celltype
sce1=subset(sce,downsample=2000)
sp1.t <- DotPlot(sce1, 
                 features =unique(genes),
                 
                 cols =Dot1c,
                 assay='RNA'  )  + coord_flip()+
  labs(x="",y="",title = "")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "T cells")+
  th
sp1.t

# NK 
Idents(sce.all)=sce.all$celltype_major
sce=subset(sce.all,idents=c('NK'))
Idents(sce)=sce$celltype
sce$celltype=factor(sce$celltype,levels=c(
  'NK_c01-CD56high', 'NK_c02-CD56low'))
genes= c(
  'GNLY','NKG7','TYROBP','KLRF1','KLRD1', 'CX3CR1','PRF1',
  'NCAM1'
)
sp1.nk <- DotPlot(sce, 
                  features =unique(genes),
                  cols =Dot1c,
                  assay='RNA'  )  + coord_flip()+
  labs(x="",y="",title = "")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "NK cells")+
  th
sp1.nk

# Myeloid 
Idents(sce.all)=sce.all$celltype_major
sce=subset(sce.all,idents=c('Mono','Platelet','DC'))
Idents(sce)=sce$celltype
sce$celltype=as.vector(sce$celltype)
sce$celltype <- factor(sce$celltype, levels = c(
  'Mono_c01-CD14-VCAN','Mono_c02-CD14-CCL3','Mono_c03-CD14-HLA-DPB1','Mono_c04-CD14-CD16',
  'DC-mDC', 'Platelet'))
Idents(sce)=sce$celltype

genes= c(
  'CD68', 'CD163', 'CD14', 
  'VCAN',
  "CCL3",
  'HLA-DPB1',
  'FCGR3A',
  'FCGR3B',
  'CD1E','CD1C', # mDC
  # platelet
  'PF4','PPBP'
  
)
sp1.m <- DotPlot(sce, 
                 features =unique(genes),
                 cols =Dot1c,
                 assay='RNA'  )  + coord_flip()+
  labs(x="",y="",title = "")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(title = "Myeloid cells")+
  th
sp1.m

layout <- c(
  area(t=2,l= 1,b=3,r=2),
  area(t=2,l= 3,b=3,r=6),
  area(t=4,l= 1,b=4,r=1),
  area(t=4,l= 3,b=5,r=6)
)

sp1_1=sp1.b+sp1.t+sp1.nk+sp1.m+
  plot_layout(design = layout) 
sp1_1