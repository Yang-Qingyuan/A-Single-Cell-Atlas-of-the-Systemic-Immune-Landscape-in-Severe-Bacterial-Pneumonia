library(Seurat)
library(ggplot2)
library(scales)
th=theme(text = element_text(
  color="black"),   
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1), 
  axis.ticks = element_line(linewidth = 0.5) , 
  axis.text = element_text(size = 10, face = "bold", color = "black"))

filep="./MiddleFile/"
outl=c('Control','Alive','Deceased')

sce.all=readRDS(paste0(filep,'4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))

## Figure 4A
colm=c(#Mono        
  '#e68b81','#f7df87','#8ccbea', '#cceac4')
library(monocle)
mycds=readRDS("pseudotime/mycds-mono_VariableFeatures.rds")

ctcol=colm
df=pData(mycds) 

df$scaled_pseudotime <- rescale(df$Pseudotime, to = c(0, 100))
df$group=df$outcome
table(df$group)

p4.1=ggplot(df, aes(x = scaled_pseudotime,color = celltype,fill = celltype)) +
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
  facet_wrap(~celltype, ncol = 1
             # ,
             # scales = "free"
  )+
  th
p4.1

## Figure 4B
Idents(sce.all)=sce.all$celltype_major
sce=subset(sce.all,idents=c('Mono'))
cell_types <- FetchData(sce, 
                        vars = c("celltype","outcome")) 
cell_types$outcome=factor(cell_types$outcome,levels = outl)

p4.2=ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x =outcome, fill =celltype), position = "fill", width = 0.75) +
  scale_fill_manual(values =colm) +
  #facet_wrap(~tissue)+
  #coord_flip()+
  labs(y ="cell proportions",x="",title = "")+
  guides(fill=guide_legend(title= ""))+
  theme_classic()+
  theme(
    strip.background = element_blank(), # 去掉分面标题的方框
    strip.text = element_text(size = 10, face = "bold"), # 设置分面标题的字体大小为10且加粗
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # 将X轴的文字垂直显示
  ) +
  th
p4.2

# Figure 4C
library(GseaVis) 
ctt=c('Mono_c01-CD14-VCAN','Mono_c02-CD14-CCL3','Mono_c03-CD14-HLA-DPB1','Mono_c04-CD14-CD16')
i=1
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0002521','GO:1903131')

p4.3=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            subPlot = 2,
            rmHt = T,
            legend.position = c(0.8,0.9),
            addPval = T,
            pvalX = 0.05,pvalY = 0.05,
            curveCol = c('#f6c8a8','#e68b81'))
p4.3


# Figure 4D-G
i=1
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0001816','GO:0007249','GO:0050729','GO:0050870')

p4.4=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            subPlot = 2,
            rmHt = T,
            legend.position = c(0.75,0.8),
            addPval = T,
            pvalX = 0.05,pvalY = 0.0,
            curveCol = c('#e68b81','#f7df87','#8ccbea', '#cceac4')
)
p4.4

i=2
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0005126','GO:0005125','GO:0002675','GO:0002696')
p4.5=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            subPlot = 2,
            rmHt = T,
            legend.position = c(0.75,0.8),
            addPval = T,
            pvalX = 0.05,pvalY = 0.0,
            curveCol = c('#e68b81','#f7df87','#8ccbea', '#cceac4'))
p4.5


i=3
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0002399','GO:0019882','GO:0002548','GO:0042110')
p4.6=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            subPlot = 2,
            rmHt = T,
            legend.position = c(0.75,0.8),
            addPval = T,
            pvalX = 0.05,pvalY = 0.0,
            curveCol = c('#e68b81','#f7df87','#8ccbea', '#cceac4')
)
p4.6

i=4
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0006397','GO:0008380','GO:0005681')
p4.7=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            subPlot = 2,
            rmHt = T,
            legend.position = c(0.75,0.8),
            addPval = T,
            pvalX = 0.05,pvalY = 0.0,
            curveCol = c('#f7df87','#8ccbea', '#cceac4')
)
p4.7