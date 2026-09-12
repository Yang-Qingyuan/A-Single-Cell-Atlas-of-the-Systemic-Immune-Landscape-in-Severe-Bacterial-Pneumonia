library(rio)
library(dplyr)
library(tidyr)
library(ggplot2)

filep <- "./MiddleFile/"

th <- theme(text = element_text(color = "black"),
            axis.title = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 12, face = "bold"),
            plot.caption = element_text(size = 8),
            axis.line = element_line(linewidth = 1),
            axis.ticks = element_line(linewidth = 0.5),
            axis.text = element_text(size = 10, face = "bold", color = "black"))

hmc <- c('#8ccbea', '#f7df87', '#e68b81')

## Figure S6C
treg <- import(paste0(filep, "FigS6C_CD4_c04_FOXP3_DotPlot_SourceData.xlsx"))

gene_order <- c(
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
  "TNFRSF18", "TNFRSF4", "IL10", "TGFB1",
  "ENTPD1", "NT5E", "LAG3", "HAVCR2"
)

treg$Gene <- factor(treg$Gene, levels = gene_order)
treg$Outcome <- factor(treg$Outcome, levels = c("Control", "Alive", "Deceased"))

ggplot(treg, aes(x = Gene, y = Outcome)) +
  geom_point(aes(size = Percent_expressed, color = Scaled_average_expression)) +
  scale_color_gradient(low = hmc[1], high = hmc[3]) +
  coord_flip() +
  labs(x = '', y = '', title = 'T_CD4_c04-FOXP3') +
  theme_classic() +
  th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Figure S6D
meta <- import(paste0(filep, 'meta-T_TCRmerged.xlsx'))
meta <- meta[meta$clonalFrequency != 0, ]
meta$clonalStatus <- ifelse(meta$clonalFrequency %in% c(1), "Nonexpanded", 'Expanded')

proportions <- meta %>%
  group_by(outcome, celltype, clonalStatus) %>%
  summarize(count = n()) %>%
  spread(key = clonalStatus, value = count, fill = 0) %>%
  mutate(expanded_proportion = Expanded / (Expanded + Nonexpanded))

proportions_wide <- proportions %>%
  filter(outcome %in% c("Alive", "Deceased")) %>%
  pivot_wider(names_from = outcome, values_from = expanded_proportion, values_fill = 0)

plotdf <- data.frame(celltype = proportions_wide[1:9, 1],
                     Alive = proportions_wide[1:9, 4],
                     Deceased = proportions_wide[10:18, 5])

colt <- c('#fcedbe', '#fab378', '#f49512', '#ffd377',
          '#6da5ad', '#d9e773', '#ccd4b5', '#889851', '#BAD8D2')
names(colt) <- c('T_CD4_c01-LEF1', 'T_CD4_c02-AQP3', 'T_CD4_c03-FOS', 'T_CD4_c04-FOXP3',
                 'T_CD8_c01-LEF1', 'T_CD8_c02-GZMK', 'T_CD8_c03-TYROBP', 'T_CD8_c04-ZNF683', 'T_CD8_c05-RORA')

ggplot(plotdf, aes(x = Alive, y = Deceased, color = celltype, label = celltype)) +
  geom_point(size = 3) +
  geom_text(hjust = 0, vjust = -1, size = 3.5, fontface = "bold") +
  scale_color_manual(values = colt) +
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