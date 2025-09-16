library(rio)
library(ggplot2)
library(dplyr)
library(tidyr)
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

meta=import(paste0(filep,'meta-T_TCRmerged.xlsx'))
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
plotdf=data.frame(celltype=proportions_wide[1:9,1],
                  Alive=proportions_wide[1:9,4],
                  Deceased=proportions_wide[10:18,5])


ggplot(plotdf, aes(x = Alive, y = Deceased, color = celltype, label = celltype)) +
  geom_point(size = 3) +
  geom_text(hjust = 0, vjust = -1,size = 3.5, fontface = "bold") +
  scale_color_manual(values = colt) + 
  labs(x = "Proportion of clonally expanded \n  cells in Alive", y = "Proportion of clonally expanded \n  cells in Deceased", 
       title = "") + 
  theme_minimal() +
  theme(legend.position = "none") +  
  th+
  coord_equal() +  
  xlim(0, 1) +  
  ylim(0, 1) +   
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgray")  
