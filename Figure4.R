library(Seurat)
library(ggplot2)
library(GseaVis)

filep="./MiddleFile/"

th=theme(text = element_text(color="black"),
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1),
  axis.ticks = element_line(linewidth = 0.5),
  axis.text = element_text(size = 10, face = "bold", color = "black"))

outl=c('Control','Alive','Deceased')
outcomeCol=c('#7ac7e2','#f9d580','#e3716e')
names(outcomeCol)=outl
colm=c('#f7df87','#8ccbea','#cceac4')

sce=readRDS(paste0(filep,"Mono-celltype.rds"))
sce$celltype=as.vector(sce$celltype)

## Figure 4A
cell_types <- FetchData(sce, vars = c("celltype","outcome"))
cell_types$outcome=factor(cell_types$outcome,levels = rev(outl))

ggplot(data = cell_types) +
  geom_bar(mapping = aes(x = outcome, fill = celltype), position = "fill", width = 0.75) +
  scale_fill_manual(values = colm) +
  labs(y = "cell proportions", x = "", title = "") +
  guides(fill = guide_legend(title = "")) +
  theme_classic() +
  theme(
    legend.position = 'top',
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  th +
  coord_flip()

## Figure 4C
ctt=c('Mono_c01-Classical','Mono_c02-Intermediate','Mono_c03-Non_Classical')
i=1
load(paste0(filep,ctt[i],"-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0031663','GO:0071347')

gseaNb(object = gseaRes,
       geneSetID = geneSetID,
       subPlot = 2,
       rmHt = T,
       legend.position = c(0.8,0.9),
       addPval = T,
       pvalX = 0.05, pvalY = 0.05,
       curveCol = c('#e68b81','#8ccbea','#cceac4'))

## Figure 4F
i=2
load(paste0(filep,ctt[i],"-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0002399','GO:0002478','GO:0050870')

gseaNb(object = gseaRes,
       geneSetID = geneSetID,
       subPlot = 2,
       rmHt = T,
       legend.position = c(0.8,0.9),
       addPval = T,
       pvalX = 0.05, pvalY = 0.05,
       curveCol = c('#e68b81','#f7df87','#8ccbea','#cceac4'))

## Figure 4H
Idents(sce)=sce$celltype
sce=subset(sce,idents=c("Mono_c03-Non_Classical"))
DotPlot(sce,
        features = c("FCGR3A","MS4A7","CX3CR1","LILRB1","LILRB2","IFITM3"),
        group.by = "outcome",
        idents = "Mono_c03-Non_Classical",
        cols = c('#8ccbea','#e68b81')) +
  coord_flip()
