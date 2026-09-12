library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rio)
library(ggpubr)
library(patchwork)
library(ggrepel)

filep="./MiddleFile/"

th=theme(text = element_text(color="black"),
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 0.5),
  axis.ticks = element_line(linewidth = 0.25),
  axis.text = element_text(size = 10, face = "bold", color = "black"))

colb=c('#f9d3e3','#cb7395','#8f476d',"#d35b7e")

sce=readRDS(paste0(filep,'B_Plasma-celltype.rds'))
sce$celltype=as.vector(sce$celltype)
sce$celltype=factor(sce$celltype,levels=c('B_c01-mature naive B','B_c02-memory B',
                                          'B_c03-active Plasma','B_c04-terminal Plasma'))
sce$outcome=factor(sce$outcome,levels = c("Control",'Alive','Deceased'))

## Figure 5A
cell_types <- FetchData(sce, vars = c("celltype","outcome"))

ggplot(data = cell_types) +
  geom_bar(mapping = aes(x = outcome, fill = celltype), position = "fill", width = 0.75) +
  scale_fill_manual(values = colb) +
  labs(y = "cell proportions", x = "", title = "") +
  theme_classic2() +
  guides(fill = guide_legend(title = "", nrow = 5)) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  th

## Figure 5C
Idents(sce)=sce$celltype
genes=c("PRDM1","IRF4","XBP1","MZB1",
        "TNFRSF17","JCHAIN",
        "IGHG1","IGHA1","IGKC")

DotPlot(sce, features = genes, group.by = "outcome",
            idents = "B_c04-terminal Plasma") +
  coord_flip() +
  labs(x='', y='', title='')

## Figure 5E
plot_data=import(paste0(filep,'Fig5E_PRDM1_activity_SourceData.xlsx'),
                 which = 'Violin_plot_data')
plot_data$outcome=factor(plot_data$Outcome,levels = c("Control","Alive","Deceased"))

ggplot(plot_data, aes(x = outcome, y = PRDM1_activity, fill = outcome)) +
  geom_violin(scale = "width", width = 0.8, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c(
    "Control" = "#7ac7e2",
    "Alive" = "#fabb6e",
    "Deceased" = "#b7282e"
  )) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(title = "", x = "", y = "PRDM1 target gene score") + th

## Figure 5F
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
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) + th

## Figure 5G
df=import(paste0(filep,'BCR_clonalDiversity-sample.xlsx'))
m1=sce@meta.data
m1=m1[!duplicated(m1$orig.ident),c('orig.ident','outcome')]
df1=merge(df,m1,by.x='Sample',by.y='orig.ident')

pb3.1=ggboxplot(data = df1,
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
  labs(title = 'Shannon', x = '', y = 'Clonal Diversity') +
  theme_classic() +
  th +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

pb3.2=ggboxplot(data = df1,
                x = 'outcome',
                y = 'Inv.Simpson',
                color = 'gray40',
                fill = NA,
                outlier.shape = NA
) +
  geom_jitter(aes(x = outcome, y = Inv.Simpson, color = outcome),
              width = 0.2, size = 2, alpha = 1) +
  scale_color_manual(values = c('Control' = '#7ac7e2',
                                'Alive' = '#f9d580',
                                'Deceased' = '#e3716e')) +
  labs(title = 'Inv.Simpson', x = '', y = 'Clonal Diversity') +
  theme_classic() +
  th +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

pb3=(pb3.1|pb3.2)
pb3

## Figure 5H
meta=import(paste0(filep,'meta-B_BCRmerged.xlsx'))
bcrcol=c("#cbe3e6","#efc661","#eda150","#f2cecc","#db7686","#932925")

cell_types <- meta
cell_types$outcome=factor(cell_types$outcome,levels=c("Control","Alive","Deceased"))

cell_types=cell_types[cell_types$cloneSize !="None",]
cell_types$size=ifelse(cell_types$clonalFrequency==1,"unique",
                       ifelse(cell_types$clonalFrequency <= 5, "2-5",
                              ifelse(cell_types$clonalFrequency<=20,"6-20",
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

ggplot(data = cell_types_agg, mapping = aes(x = "", y = prop, fill = size)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = bcrcol) +
  labs(y = "", x = "", title = "") +
  theme_classic() +
  th +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 10, face = 'bold'),
        strip.background = element_blank()) +
  geom_label_repel(aes(label = scales::percent(prop, accuracy = 0.01)),
                   position = position_stack(vjust = 0.5),
                   box.padding = 0.5,
                   size = 3,
                   show.legend = FALSE) +
  facet_wrap(~ outcome, nrow = 1)

## Figure 5I
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

plotdf=data.frame(celltype=proportions_wide[1:4,1],
                  Alive=proportions_wide[1:4,4],
                  Deceased=proportions_wide[5:8,5])

ggplot(plotdf, aes(x = Alive, y = Deceased, color = celltype, label = celltype)) +
  geom_point(size = 3) +
  geom_text(hjust = 0, vjust = -1, size = 3.5, fontface = "bold") +
  scale_color_manual(values = colb) +
  labs(x = "Proportion of clonally expanded \n  cells in Alive",
       y = "Proportion of clonally expanded \n  cells in Deceased",
       title = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  th +
  coord_equal() +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgray")