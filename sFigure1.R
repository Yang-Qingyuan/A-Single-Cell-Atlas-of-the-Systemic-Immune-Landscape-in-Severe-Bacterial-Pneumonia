library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(rio)

filep <- "./MiddleFile/"

th <- theme(text = element_text(family = "Arial"),
            axis.title = element_text(family = "Arial", size = 10, face = "bold"),
            plot.title = element_text(family = "Arial", size = 12, face = "bold"),
            plot.caption = element_text(family = "Arial", size = 8),
            axis.line = element_line(linewidth = 1),
            axis.ticks = element_line(linewidth = 0.5),
            axis.text = element_text(family = "Arial", size = 10, face = "bold", color = "black"))

outcomeCol <- c('#7ac7e2', '#f9d580', '#e3716e')
outl <- c('Control', 'Alive', 'Deceased')
names(outcomeCol) <- outl
sexcol <- c("#f8cbe0", "#97ceff")
outcomeCol2 <- c("#7ac7e2", '#fab378')
group_levels <- c("Control", "Severe pneumonia")

sce.all <- readRDS(paste0(filep, 'R1_4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))
meta <- sce.all@meta.data

## Figure S1A
a <- meta[!duplicated(meta$orig.ident), c("age", "outcome")]
a$age <- as.numeric(a$age)
a$outcome <- factor(a$outcome, levels = outl)
df <- a

s_p1 <- ggplot(df, aes(x = outcome, y = age, color = outcome)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.18, size = 2, alpha = 0.8) +
  scale_color_manual(values = outcomeCol) +
  scale_x_discrete(limits = outl) +
  theme_classic(base_size = 14) +
  labs(title = "scRNA-seq",
       x = NULL, y = "Age") +
  theme(legend.position = "none") + th


b <- import(paste0(filep, 'flow_cytometry_cohort_clinical_info.xlsx'))
b$age <- as.numeric(b$Age)
b$outcome <- factor(b$Outcome, levels = outl)
df <- b

f_p1 <- ggplot(df, aes(x = outcome, y = age, color = outcome)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.18, size = 2, alpha = 0.8) +
  scale_color_manual(values = outcomeCol) +
  scale_x_discrete(limits = outl) +
  theme_classic(base_size = 14) +
  labs(title = 'scRNA-seq',
       x = NULL, y = "Age") +
  theme(legend.position = "none") + th


## Figure S1B
a <- meta[!duplicated(meta$orig.ident), ]
count_data <- a %>%
  group_by(outcome, sex) %>%
  summarise(count = n()) %>%
  ungroup()
count_data$outcome <- factor(count_data$outcome, levels = outl)

sp1.2 <- ggplot(count_data, aes(x = outcome, y = count, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = sexcol) +
  labs(title = 'scRNA-seq',
       x = "", y = "Number  of Donors", fill = "Sex") +
  theme_classic() +
  th +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title = NULL))

count_data <- b %>%
  group_by(outcome, Sex) %>%
  summarise(count = n()) %>%
  ungroup()
count_data$outcome <- factor(count_data$outcome, levels = outl)

f_p1.2 <- ggplot(count_data, aes(x = outcome, y = count, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = sexcol) +
  labs(title = 'Flow cytometry',
       x = "", y = "Number  of Donors", fill = "Sex") +
  theme_classic() +
  th +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title = NULL))


## Figure S1C
sp1.3 <- ggplot(meta, aes(x = outcome, y = nCount_RNA, fill = outcome)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Number of Counts") +
  theme_classic() +
  theme(legend.position = "none") + th
sp1.3

## Figure S1D
sp1.4 <- ggplot(meta1, aes(x = outcome, y = nFeature_RNA, fill = outcome)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Number of Genes") +
  theme_classic() +
  theme(legend.position = "none") + th
sp1.4

## Figure S1E
sp1.5 <- ggplot(meta1, aes(x = outcome, y = percent_mito, fill = outcome)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Percentage of \n Mitochondrial genes") +
  theme_classic() +
  theme(legend.position = "none") + th
sp1.5


## Figure S1F-I
Dot1c <- c("#F8D941", "#C73866")

Idents(sce.all) <- sce.all$celltype_major
sce <- subset(sce.all, idents = c("B", "Plasma"))
Idents(sce) <- sce$celltype

sce$celltype <- factor(sce$celltype, levels = c('B_c01-mature naive B', 'B_c02-memory B',
                                                'B_c03-active Plasma', 'B_c04-terminal Plasma'))
genes <- c(
  'CD19', 'CD79A', 'CD79B', 'MS4A1', 'CD20',
  'IGHD', 'IGHM', 'CCR7', 'TCL1A', 'FAM129C',
  'EBI3', 'FCRL4', 'DUSP4', 'GPR183',
  'PRDM1', 'XBP1', 'CCR10', 'TXNDC5',
  'LGALS3', 'CCL2', 'CCL5', 'ANXA1'
)

sp1.b <- DotPlot(sce,
                 features = unique(genes),
                 cols = Dot1c,
                 assay = 'RNA') + coord_flip() +
  labs(x = "", y = "", title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "B cells") +
  th
sp1.b

Idents(sce.all) <- sce.all$celltype_major
sce <- subset(sce.all, idents = c('Mono', 'Platelet', 'DC'))
Idents(sce) <- sce$celltype

sce$celltype <- as.vector(sce$celltype)
sce$celltype <- factor(sce$celltype, levels = c('Mono_c01-Classical', 'Mono_c02-Intermediate', 'Mono_c03-Non_Classical',
                                                'DC-mDC', 'Platelet'))
Idents(sce) <- sce$celltype

genes <- c(
  'CD68', 'CD163', 'CD14',
  'VCAN',
  'HLA-DPB1', 'CD74',
  'HLA-DRA', 'HLA-DRB1',
  'FCGR3A',
  'FCGR3B',
  'CX3CR1',
  'CD1E', 'CD1C',
  'PF4', 'PPBP',
  'CSF3R', 'MPO', 'ELANE'
)

sp1.m <- DotPlot(sce,
                 features = unique(genes),
                 cols = Dot1c,
                 assay = 'RNA') + coord_flip() +
  labs(x = "", y = "", title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Myeloid cells") +
  th
sp1.m

Idents(sce.all) <- sce.all$celltype_major
sce <- subset(sce.all, idents = c('CD4+T', 'CD8+T'))
sce$celltype <- as.vector(sce$celltype)
sce$celltype <- factor(sce$celltype, levels = c(
  'T_CD4_c01-LEF1', 'T_CD4_c02-AQP3', 'T_CD4_c03-FOS', 'T_CD4_c04-FOXP3',
  'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK', 'T_CD8_c03-TYROBP', 'T_CD8_c04-ZNF683', 'T_CD8_c05-RORA'
))

genes <- c(
  'CD3D', 'CD3E',
  'CD4',
  'LEF1', 'AQP3', 'FOS', 'FOXP3',
  'CD8A', 'CD8B',
  'GZMK', 'TYROBP', 'ZNF683', 'RORA'
)

Idents(sce) <- sce$celltype
sce1 <- subset(sce, downsample = 2000)
sp1.t <- DotPlot(sce1,
                 features = unique(genes),
                 cols = Dot1c,
                 assay = 'RNA') + coord_flip() +
  labs(x = "", y = "", title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "T cells") +
  th
sp1.t

Idents(sce.all) <- sce.all$celltype_major
sce <- subset(sce.all, idents = c('NK'))
Idents(sce) <- sce$celltype
sce$celltype <- factor(sce$celltype, levels = c(
  'NK_c01-CD56high', 'NK_c02-CD56low'))
genes <- c(
  'GNLY', 'NKG7', 'TYROBP', 'KLRF1', 'KLRD1', 'CX3CR1', 'PRF1',
  'NCAM1'
)
sp1.nk <- DotPlot(sce,
                  features = unique(genes),
                  cols = Dot1c,
                  assay = 'RNA') + coord_flip() +
  labs(x = "", y = "", title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "NK cells") +
  th
sp1.nk

## Figure S1J
ctl <- c('B_c01-mature naive B', 'B_c02-memory B',
         'B_c03-active Plasma', 'B_c04-terminal Plasma',
         'T_CD4_c01-LEF1', 'T_CD4_c02-AQP3', 'T_CD4_c03-FOS', 'T_CD4_c04-FOXP3',
         'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK', 'T_CD8_c03-TYROBP', 'T_CD8_c04-ZNF683', 'T_CD8_c05-RORA',
         'NK_c01-CD56high', 'NK_c02-CD56low',
         'DC-mDC',
         'Mono_c01-Classical', 'Mono_c02-Intermediate', 'Mono_c03-Non_Classical',
         'Platelet'
)
ctc <- c('#f9d3e3', '#cb7395',
         '#8f476d', "#d35b7e",
         '#fcedbe', '#fab378', '#f49512', '#ffd377',
         '#6da5ad', '#d9e773',
         '#ccd4b5', '#889851',
         '#BAD8D2',
         "#c9d4f7", '#caadd8',
         '#b1ebfe',
         '#d7cdbe', '#6B4C3D', '#cd9e80',
         '#b9181a')
names(ctc) <- ctl

ctl_big <- c('B', 'Plasma', 'CD4+T', 'CD8+T', 'NK', 'DC', 'Mono', 'Platelet')
ctc_big <- c('#f091a0', "#d35b7e", '#fab378', '#80b5b8', '#caadd8', '#80b1d2', '#6B4C3D', '#b9181a')
sce.all$celltype_major <- factor(sce.all$celltype_major, levels = ctl_big)
sce.all$celltype <- as.vector(sce.all$celltype)
sce.all$celltype <- factor(sce.all$celltype, levels = ctl)

Idents(sce.all) <- sce.all$outcome
con <- subset(sce.all, idents = c('Control'))
ali <- subset(sce.all, idents = c('Alive'))
dec <- subset(sce.all, idents = c('Deceased'))

DimPlot(con, group.by = "celltype", reduction = "umap_harmony",
        raster = TRUE,
        cols = ctc) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '') + th +
  theme(legend.position = 'bottom', legend.justification = 'center') +
  guides(color = guide_legend(ncol = 2))

DimPlot(ali, group.by = "celltype", reduction = "umap_harmony",
        raster = TRUE,
        cols = ctc) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '') + th +
  theme(legend.position = 'bottom', legend.justification = 'center') +
  guides(color = guide_legend(ncol = 2))

DimPlot(dec, group.by = "celltype", reduction = "umap_harmony",
        raster = TRUE,
        cols = ctc) +
  labs(x = 'UMAP1', y = 'UMAP2', title = '') + th +
  theme(legend.position = 'bottom', legend.justification = 'center') +
  guides(color = guide_legend(ncol = 2))
