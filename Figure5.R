library(ggplot2)
library(Seurat)
filep="./MiddleFile/"

th=theme(text = element_text(
  color="black"), 
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 0.5), 
  axis.ticks = element_line(linewidth =0.25) ,
  axis.text = element_text(size = 10, face = "bold", color = "black"))


colb=c('#f9d3e3','#cb7395','#ffadbb','#8f476d',"#d35b7e")
ctcol=colb
groupcol=c('#7ac7e2','#fabb6e','#b7282e')
colbd=c("#F8D941","#C73866")

## Figure 5A-B
library(monocle)
library(scales)
mycds=readRDS("pseudotime/mycds-B.rds")
df=pData(mycds) 

df$scaled_pseudotime <- rescale(df$Pseudotime, to = c(0, 100))
df$group=df$outcome
df$group=factor(df$group,levels=c("Control","Alive",'Deceased'))

ggplot(df, aes(x = scaled_pseudotime,color = celltype,fill = celltype)) +
  geom_density(alpha = 0.5, linewidth = 1) +  
  theme_classic() +
  scale_color_manual(values = ctcol) +
  scale_fill_manual(values = ctcol) +
  labs(x = "Pseudotime (scaled)", y = "") +  
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none", 
    strip.text = element_blank(), 
    strip.background = element_blank(),  
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),  
    axis.text.x = element_text(angle = 0, hjust = 0.5) 
  ) +
  facet_wrap(~celltype, ncol = 1,
             scales = "free")+
  th

df1=df[df$celltype=='B_c05-terminal Plasma',]
ggplot(df1, aes(x = scaled_pseudotime,color = group,fill = group)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  theme_classic() +
  scale_color_manual(values = groupcol) +
  scale_fill_manual(values = groupcol) +
  labs(x = "Pseudotime (scaled)", y = "Density") + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.title = element_blank(),  
    strip.text = element_blank(), 
    strip.background = element_blank(), 
    # axis.text.y = element_blank(),  
    # axis.ticks.y = element_blank(),  
    axis.text.x = element_text(angle = 0, hjust = 0.5),  
    legend.position = "right" 
  ) +
  th

## Figure 5C
sce=readRDS(paste0(filep,'B_Plasma-celltype.rds'))
sce$celltype=as.vector(sce$celltype)
sce$celltype=factor(sce$celltype,levels=c('B_c01-mature naive B','B_c02-memory B',
                                          'B_c03-Plasmablast',
                                          'B_c04-active Plasma','B_c05-terminal Plasma'))
sce$outcome=factor(sce$outcome,levels = c("Control",'Alive','Deceased'))
cell_types <- FetchData(sce, 
                        vars = c("celltype","outcome")) 
ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x =outcome, fill =celltype), position = "fill", width = 0.75) +
  scale_fill_manual(values =colb) +
  labs(y ="cell proportions",x="",title = "")+
  guides(fill=guide_legend(title= "Cell type"))+
  theme_classic()+
  guides(fill = guide_legend(title = "", nrow = 5)) + 
  theme_classic()+
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  th

## Figure 5D
Idents(sce)=sce$celltype
sce5=subset(sce,idents=c('B_c05-terminal Plasma'))
genes=c(
  "PRDM1", "IRF4", "XBP1", "MZB1", 
  "TNFRSF17", "JCHAIN", 
  "IGHG1", "IGHA1", "IGKC")

pb1=DotPlot(sce5, features = genes, group.by = "outcome")+
  coord_flip() +
  labs(x='',y='',title='')
data1 <-pb1$data
combined_data <-data1
pct_exp_range <- range(combined_data$pct.exp, na.rm = TRUE)
avg_exp_scaled_range <- range(combined_data$avg.exp.scaled, na.rm = TRUE)

ggplot(data1, aes(x =id , y = features.plot, size = pct.exp,
                          fill = avg.exp.scaled)) +
  geom_point(shape = 21, stroke = 0.5, color = "black") +
  scale_fill_gradient2(low = "#79c4ff", mid="#fffb7d", high ="#e86898", 
                       name = "Mean expression\nin group (mRNA)", 
                       limits = avg_exp_scaled_range) +
  scale_size_continuous(name = "Fraction of cells\nin group (%)", 
                        range = c(2, 8), 
                        limits = pct_exp_range) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0), 
    plot.title.position = "plot", 
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "", y = "", title ='' )+th

## Figure 5E-F
library(circlize)
library(rio)
colors <- colorRampPalette(c(
  '#7ac7e2','#f9d580','#e3716e'))(100)
values <- seq(0, 303, length.out = 101)[-101]
col_fun <- colorRamp2(values, colors)
col_fun

   ## upregulate
pathn=c('B cell receptor signaling pathway','adaptive immune response','plasma membrane invagination','complement activation, classical pathway','humoral immune response mediated by circulating immunoglobulin','immunoglobulin mediated immune response')
a=import(paste0(filep,"Deceased-TerminialPlasma-gseGO.xlsx"))
colnames(a)
a=a[a$Description %in% pathn,]
df1=a[match(pathn,a$Description),]
           
ggplot(df1, aes(x = reorder(Description, NES), y = NES, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = col_fun(values)) +
  theme_classic() +
  labs(x = "", y = "NES", fill = "-log10(p.adjust)", title = '') +
  coord_flip() +  
  theme(
    legend.title = element_text(size = 9,face = 'bold'),  
    legend.position = 'right'
  )+th

  ## downregulate
pathn=c('positive regulation of defense response','positive regulation of response to external stimulus',
'positive regulation of inflammatory response','positive regulation of interleukin-1 beta production')
a=import(paste0(filep,"Deceased-TerminialPlasma-gseGO.xlsx"))
colnames(a)
a=a[a$Description %in% pathn,]
df1=a[match(pathn,a$Description),]

ggplot(df1, aes(x = reorder(Description, NES), y = NES, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = col_fun(values)) +
  theme_classic() +
  labs(x = "", y = "NES", fill = "-log10(p.adjust)", title = '') +
  coord_flip() +  # Flip coordinates
  theme(
    #axis.text.y = element_markdown(size = 9),  # Set Y-axis text size to 9
    legend.title = element_text(size = 9,face = 'bold'),  # Set legend title font size to 9
    legend.position = 'right'
  )+th

## Figure 5G
library(dplyr)
meta=import(paste0(filep,'meta-BCRisotype.xlsx'))
bcricol=c("#b5dcf3","#71b8ec",
                   "#c1e2db","#82c4b8",
                   "#f9dbb5","#f1b56d",'#f9bdb5','#f67c6f')

meta$BCR_isotype=factor(meta$BCR_isotype,levels = 
                          c('IGHM','IGHD','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4'))
meta$outcome=factor(meta$outcome,levels = c('Control','Alive','Deceased'))

meta_agg <- meta %>%
  group_by(outcome, BCR_isotype) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

ggplot(data = meta_agg, mapping = aes(x = outcome, y = percentage, fill = BCR_isotype)) +
  geom_bar(stat = "identity", 
           position = "fill", 
           width = 0.7) +
  scale_fill_manual(values = bcricol) +
  labs(y = "BCR isotype (%)", x = "", title = "", fill = "") +
  theme_classic()+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold")
  )+th

## Figure 5H
library(ggpubr)
library(patchwork)
combined=readRDS(paste0(filep,'combinedBCR-AddMeta_Sample.rds'))
df=as.data.frame(clonalDiversity(combined, 
                                 cloneCall = "strict", 
                                 group.by = "Sample",
                                 exportTable=TRUE,
                                 n.boots = 100))
export(df,paste0(filep,'BCR_clonalDiversity-sample.xlsx'))

df=import(paste0(filep,'BCR_clonalDiversity-sample.xlsx')) 
m1=import(paste0(filep,'Control_Alive_Deceased-meta.xlsx'))
df1=merge(df,m1,by.x='Sample',by.y='orig.ident')

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
                     label = 'p.signif') +
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
                     label = 'p.signif') +
  labs(title = '', x = '', y = 'Inv.Simpson index') +
  theme_classic()+
  th

## Figure 5I
meta=import(paste0(filep,'meta-B_BCRmerged.xlsx'))
bcrcol=c("#cbe3e6",
         "#efc661","#eda150","#f2cecc","#db7686","#932925")

cell_types <- meta
cell_types$outcome=factor(cell_types$outcome,levels=c("Control","Alive","Deceased"))

cell_types=cell_types[cell_types$cloneSize !="None",]
cell_types$size=ifelse(cell_types$clonalFrequency==1,"unique",
                       ifelse(cell_types$clonalFrequency <= 5, "2-5",ifelse(cell_types$clonalFrequency<=20,"6-20",
                                                                                                                           ifelse(cell_types$clonalFrequency<=50,"21-50",
                                                                                                                                  ifelse(cell_types$clonalFrequency<=100,"51-100", ">100")))))

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

library(ggrepel)
ggplot(data = cell_types_agg, mapping = aes(x = "", y = prop, fill = size)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = bcrcol) +
  labs(y = "", x = "", title = "") +
  theme_classic() +
  th+
  theme(legend.position = "right", 
        legend.text = element_text(size = 8), 
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),   
        axis.ticks = element_blank(),    
        axis.line = element_blank(),   
        strip.text = element_text(size = 10,face = 'bold'), 
        strip.background = element_blank()) +  
  geom_label_repel(aes(label = scales::percent(prop, accuracy = 0.01)),  
                   position = position_stack(vjust = 0.5),             
                   box.padding = 0.5,                                 
                   size = 3,
                   show.legend = FALSE) +                             
  facet_wrap(~ outcome, nrow = 1)  

## Figure 5J
library(tidyr)
meta=import(paste0(filep,'meta-B_BCRmerged.xlsx'))
meta=meta[meta$clonalFrequency!=0,]
meta$clonalStatus=ifelse(meta$clonalFrequency %in% c(1),"Nonexpanded",'Expanded')

proportions <- meta %>%
  group_by(outcome, celltype, clonalStatus) %>%
  summarize(count = n()) %>%
  spread(key = clonalStatus, value = count, fill = 0) %>%
  mutate(expanded_proportion = Expanded / (Expanded + Nonexpanded))

proportions_wide <- proportions %>%
  filter(outcome %in% c("Alive", "Deceased")) %>%
  pivot_wider(names_from = outcome, values_from = expanded_proportion, values_fill = 0)

plotdf=data.frame(celltype=proportions_wide[1:5,1],
                  Alive=proportions_wide[1:5,4],
                  Deceased=proportions_wide[6:10,5])

ggplot(plotdf, aes(x = Alive, y = Deceased, color = celltype, label = celltype)) +
  geom_point(size = 3) +
  geom_text(hjust = 0, vjust = -1,size = 3.5, fontface = "bold") +  
  scale_color_manual(values = colb) +  
  labs(x = "Proportion of clonally expanded \n  cells in Alive", y = "Proportion of clonally expanded \n  cells in Deceased", 
       title = "") +
  theme_minimal() +
  theme(legend.position = "none") +  
  th+
  coord_equal() + 
  xlim(0, 1) +  
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgray") 
