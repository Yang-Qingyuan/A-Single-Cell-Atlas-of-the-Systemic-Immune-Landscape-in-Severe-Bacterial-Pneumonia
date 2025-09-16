library(rio)
library(ggplot2)
filep="./MiddleFile/"

th=theme(text = element_text(
  color="black"),  
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1),  
  axis.ticks = element_line(linewidth = 0.5) ,
  axis.text = element_text(size = 10, face = "bold", color = "black"))

## sFigure 5A
meta=import(paste0(filep,'meta-BCRisotype.xlsx'))
bcricol=c("#b5dcf3","#71b8ec",
                   "#c1e2db","#82c4b8",
                   "#f9dbb5","#f1b56d",'#f9bdb5','#f67c6f')
meta$BCR_isotype=factor(meta$BCR_isotype,levels = 
                          c('IGHM','IGHD','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4'))
meta$outcome=factor(meta$outcome,levels = c('Control','Alive','Deceased'))

meta_agg <- meta %>%
  group_by(outcome, BCR_isotype,celltype) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

ggplot(data = meta_agg, mapping = aes(x = celltype, y = percentage, fill = BCR_isotype)) +
  geom_bar(stat = "identity", 
           position = "fill", 
           width = 0.7) +
  scale_fill_manual(values = bcricol) +
  labs(y =  "BCR isotype (%)", x = "", title = "", fill = "") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),   
    strip.text = element_text(face = "bold")  
  )+th

## sFigure 5B
meta=import(paste0(filep,'meta-B_BCRmerged.xlsx'))
meta=meta[meta$clonalFrequency!=0,]
bcrcol1=c("#cbe3e6","#efc661")

meta$clonestatus=ifelse(meta$clonalFrequency==1,"No clonal","Clonal")
meta$clonestatus=factor(meta$clonestatus,levels = 
                          c("No clonal","Clonal"))
meta$outcome=factor(meta$outcome,levels = c('Control','Alive','Deceased'))

meta_agg <- meta %>%
  group_by(outcome, celltype, clonestatus) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

ggplot(data = meta_agg, mapping = aes(x = celltype, y = percentage, fill = clonestatus)) +
  geom_bar(stat = "identity", 
           position = "fill", 
           width = 0.7) +
  facet_wrap(~ outcome, 
             nrow = 1) +
  scale_fill_manual(values = bcrcol1) +
  labs(y = "Distribution of clone status (%)", x = "", title = "") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold") 
  )+th
