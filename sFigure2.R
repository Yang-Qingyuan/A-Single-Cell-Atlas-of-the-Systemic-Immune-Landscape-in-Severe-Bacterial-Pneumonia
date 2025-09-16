library(Seurat)
library(rio)
library(tidyverse)
library(ggpubr)

th=theme(text = element_text(#family = "Arial",
  color="black"), 
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 0.5),  
  axis.ticks = element_line(linewidth =0.25) , 
  axis.text = element_text(size = 10, face = "bold", color = "black"))

filep="./MiddleFile/"

## sFigure 2A
df=import(paste0(filep,'meta-total.xlsx'))
colnames(df)
table(df$celltype_major)
b=df[!duplicated(df$orig.ident), c('orig.ident','outcome')]
result <- df %>%
  group_by(orig.ident, celltype_major) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
data=merge(result,b,by="orig.ident")

plot_boxplot <- function(data, ct) {
  cell_data <- data[data$celltype_major == ct, ]
  max_y <- max(cell_data$Proportion, na.rm = TRUE)
  
  p <- ggplot(cell_data, aes(x = outcome, y = Proportion, fill = outcome)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2, size = 1, alpha = 0.7) + 
    theme_classic() +  
    theme(
      axis.title.x = element_blank(),  
      plot.title = element_text(hjust = 0.5) 
    ) + 
    labs(title = ct, y = "Proportion", x = '') +  
    scale_fill_manual(values = c("Control" = '#7ac7e2', "Alive" = '#f9d580', "Deceased" = '#e3716e')) +  
    geom_signif(
      comparisons = list(
        c("Control", "Alive"),
        c("Alive", "Deceased"),
        c("Control", "Deceased")
      ),
      map_signif_level = TRUE,
      y_position = c(max_y * 1.1, max_y * 1.1, max_y * 1.2)  
    ) + 
    th
}

data$outcome=factor(data$outcome,levels = c('Control','Alive','Deceased'))
ct <- c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono','Platelet')
plots <- lapply(ct, function(cyt) plot_boxplot(data, cyt))
p1 <- wrap_plots(plots, ncol = 4)
p1