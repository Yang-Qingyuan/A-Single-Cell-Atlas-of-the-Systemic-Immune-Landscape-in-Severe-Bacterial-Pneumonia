library(Seurat)
library(ggplot2)
library(dplyr)
th=theme(text = element_text(
  color="black"),  
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1.5),  
  axis.ticks = element_line(linewidth = 1), 
  axis.text = element_text(size = 10, face = "bold", color = "black"))

colt=c( #CD4+T
  '#fcedbe','#fab378','#f49512','#ffd377',
           #CD8+T
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
sce <- readRDS(paste0(filep,'T-celltype.rds'))
sce$outcome=factor(sce$outcome,levels=c('Control','Alive','Deceased'))
cytotoxic_genes <- c(
  "PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH",
  "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7", "CCL4", "CCL3")
naive_genes <- c(
  "IL7R", "CCR7", "SELL", "FOXO1", "KLF2", "KLF3", "LEF1", "TCF7", "ACTN1", "FOXP1")
regulatory_genes <- c(
  "IL2RA", "IL1R2", "FOXP3", "LAYN", "TNFRSF9", "FANK1",
  "RTKN2", "IL1R1", "CUL9", "IKZF2")
gene_modules <- list(
  Cytotoxic = cytotoxic_genes,
  Naive = naive_genes,
  Regulatory = regulatory_genes
)
scRNAsub=AddModuleScore(sce,features = gene_modules,name = names(gene_modules))
meta <- scRNAsub@meta.data
scores_df <- meta[, c("celltype", "Cytotoxic1", "Naive2", "Regulatory3")]

data <- scores_df %>%
  group_by(celltype) %>%
  summarise(
    Cytotoxic1 = mean(Cytotoxic1, na.rm = TRUE),
    Naive2 = mean(Naive2, na.rm = TRUE),
    Regulatory3 = mean(Regulatory3, na.rm = TRUE)
  ) %>%
  as.data.frame()

colnames(data)=c("T_celltype", "Cytotoxic", "Naive", "Regulatory")

angles <- (c(0, 120, 240) * pi / 180) 

data$x <- with(data, Cytotoxic * cos(angles[1]) + Naive * cos(angles[2]) + Regulatory * cos(angles[3]))
data$y <- with(data, Cytotoxic * sin(angles[1]) + Naive * sin(angles[2]) + Regulatory * sin(angles[3]))

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
library(monocle)
mycds=readRDS("pseudotime/mycds-CD4T.rds")
ctcol=c( #CD4+T
  '#fcedbe','#fab378','#f49512','#ffd377')
df=pData(mycds)   
df$scaled_pseudotime <- rescale(df$Pseudotime, to = c(0, 100))

ggplot(df, aes(x = scaled_pseudotime,color = celltype,fill = celltype)) +
  geom_density(alpha = 0.5, linewidth = 1) +  
  theme_classic() +
  scale_color_manual(values = ctcol) +
  scale_fill_manual(values = ctcol) +
  labs(x = "Pseudotime (scaled)", y = "") + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none", 
    strip.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.x = element_text(angle = 0, hjust = 0.5) 
  ) +
  facet_wrap(~celltype, ncol = 1)+
  th

mycds=readRDS("pseudotime/mycds-CD8T.rds")
ctcol=c(#CD8+T
  '#6da5ad','#d9e773',
           '#ccd4b5','#889851', 
           '#BAD8D2')
df=pData(mycds)
df$scaled_pseudotime <- rescale(df$Pseudotime, to = c(0, 100))

ggplot(df, aes(x = scaled_pseudotime,color = celltype,fill = celltype)) +
  geom_density(alpha = 0.5, linewidth = 1) + 
  theme_classic() +
  scale_color_manual(values = ctcol) +
  scale_fill_manual(values = ctcol) +
  labs(x = "Pseudotime (scaled)", y = "") +  
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none", 
    strip.background = element_blank(), 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    axis.text.x = element_text(angle = 0, hjust = 0.5)  
  ) +
  facet_wrap(~celltype, ncol = 1)+
  th

## Figure 6C
meta=sce@meta.data
meta$celltype_major=as.vector(meta$celltype_major)
meta$celltype=as.vector(meta$celltype)
meta1=meta[!duplicated(meta$orig.ident),]
df=meta
grouped_big <- df %>%
  group_by(orig.ident, celltype_major, celltype) %>%
  summarise(count = n()) %>%
  ungroup()

total_counts_big <- df %>%
  group_by(orig.ident, celltype_major) %>%
  summarise(total_count_big = n()) %>%
  ungroup()

merged_big <- grouped_big %>%
  left_join(total_counts_big, by = c("orig.ident", "celltype_major"))

merged_big <- merged_big %>%
  mutate(proportion = count / total_count_big)

merged_big1=merged_big[,c('orig.ident','celltype','proportion')]
m2=merge(merged_big1,meta1[,c('orig.ident','outcome')])

ggboxplot(data = m2, 
          x = 'outcome', 
          y = 'proportion' ,
          color = 'gray40',   
          fill = NA,         
          outlier.shape = NA  
) +
  geom_jitter(aes(x = outcome, y = proportion, color = outcome),
              width = 0.2, size = 2, alpha = 0.5) +
  scale_color_manual(values = c('Control' = '#7ac7e2', 
                                'Alive' = '#f9d580', 
                                'Deceased' = '#e3716e')) +
  stat_compare_means(comparisons = list(c('Control','Alive'),
                                        c('Alive','Deceased'),
                                        c('Control','Deceased')),
                     label = 'p.signif') +
  labs(title = '', x = '', y = 'proportion') +
  theme_classic() +
  th +
  facet_wrap(~celltype, scales = 'free', ncol = 5) +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10), 
    axis.text.x = element_blank()  
  )

## Figure 6D-E
#PAGA analysis in Python

## Figure 6F
library(pheatmap)
meta=import(paste0(filep,'CD4_c03-score_meta.xlsx'))
meta1=meta[meta$celltype %in% c('T_CD4_c03-FOS'),]

selected_columns <-c("TCR signaling","Costimulatory molecules",
                     "Effector-primed signature",
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
hm=hm[,c(2,1,3)]
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = 'T_CD4_c03-FOS')

## Figure 6G
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
hm=hm[,c(2,1,3)]
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
hm=hm[,c(2,1,3)]
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
hm=hm[,c(2,1,3)]
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
hm=hm[,c(2,1,3)]
pheatmap(hm,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(hmc)(100),
         angle_col ='90',
         scale='row',
         main = cn)


## Figure 6H
Idents(sce)=sce$celltype
cd4_4=subset(sce,idents=c('T_CD4_c04-FOXP3'))

scorelist <- list(
  Immune_suppression= c("FOXP3", "IL2RA","IKZF2")
)
t=AddModuleScore(cd4_4,features = scorelist, search = TRUE)
colnames(t@meta.data)[(length(colnames(t@meta.data))-(length(scorelist)-1)):length(colnames(t@meta.data))]=names(scorelist)
meta=t@meta.data
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
  stat_compare_means(comparison=list(c('Control','Alive'),c('Alive','Deceased'),c('Control','Deceased')),
                     label='p.signif')+
 scale_fill_manual(values = outcomeCol) +
  labs(x = "", y = "Immune suppression",title = 'T_CD4_c04-FOXP3') +
  theme_classic() +
  theme(legend.position = "none")+th

## Figure 6I
df=rio::import(paste0(filep,'TCR_clonalDiversity-sample.xlsx'))
meta=sce@meta.data
meta1=meta[!duplicated(meta$orig.ident) ,c('orig.ident','outcome')]
df1=merge(df[,c('Shannon','Inv.Simpson','Sample')],meta1,by.x='Sample',by.y='orig.ident')

ggboxplot(data = df1, 
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
  stat_compare_means(comparisons = list(c('Control','Alive'),
                                        c('Alive','Deceased'),
                                        c('Control','Deceased')),
                     label = 'p.format'# 'p.signif'
  ) +
  labs(title = '', x = '', y = 'Shannon index') +
  theme_classic()+
  th

ggboxplot(data = df1, 
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
  stat_compare_means(comparisons = list(c('Control','Alive'),
                                        c('Alive','Deceased'),
                                        c('Control','Deceased')),
                     label = 'p.format'# 'p.signif'
  ) +
  labs(title = '', x = '', y = 'Inv.Simpson index') +
  theme_classic()+
  th

### Figure 6J
library(ggrepel)
meta=import(paste0(filep,'meta-T_TCRmerged.xlsx'))
tcrcol=c("#cbe3e6",
                  '#889851','#d9e773','#fcedbe','#ffd377','#fab378')
tcrcol1=c("#cbe3e6",'#d9e773','#fab378')

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
cell_types1=cell_types
cell_types_agg <- cell_types1 %>%
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

