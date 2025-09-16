library(Seurat)
library(ggplot2)
library(GseaVis) 
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
ctt=c('Mono_c01-CD14-VCAN','Mono_c02-CD14-CCL3','Mono_c03-CD14-HLA-DPB1','Mono_c04-CD14-CD16')

sce.all=readRDS(paste0(filep,'4.publicHC_PBMC-sce.all_qc_cluster-join-celltype.rds'))

## sFigure 4A
i=2
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0002521')
p1.1=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            #subPlot = 2,
            rmHt = T,
            termWidth = 35,
            #legend.position = c(0.8,-0.1),
            addPval = T,
            #pvalX = 0.05,pvalY = 0.05,
            curveCol = c('#f6c8a8','#e68b81'))


## sFigure 4B
i=1
load(paste0(filep,ctt[i],'-',"Deceased-gseGO.RData"))
gseaRes=GO_kk
geneSetID = c('GO:0008219','GO:0006915','GO:0012501')

p1.2=gseaNb(object = gseaRes,
            geneSetID = geneSetID,
            #subPlot = 2,
            rmHt = T,
            termWidth = 35,
            legend.position = c(0.8,0.8),
            addPval = T,
            pvalX = 0.05,pvalY = 0.05,
            curveCol = c('#c5dff4','#80c1c4','#b2d3a4'))


## sFigure 4C
library(rio)
library(pheatmap)
Idents(sce.all)=sce.all$celltype_major
mono=subset(sce.all,idents=c('Mono'))
table(mono$celltype)
mono$celltype=as.vector(mono$celltype)

scorelist=list(
  Pathogen_recognition_and_clearance=c("TLR2", "TLR4", "NOD2", "CLEC7A"),
  Pro_inflammatory=c('IL1B','IL6','TNF',#'MM9',
                     'TNFRSF12A','LPAR3',"CCR5",
                     'CCL2','CCL3','CCL4','CXCL2','CXCL9','CXCL10','SOCS3'),
  Antigen_presentation=c( "CD74", "CIITA", "CD86", "CD80"),
  
  MHC_II=c('HLA-DRA','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DPB1','HLA-DMB')
)
mono=AddModuleScore(mono,features = scorelist, search = TRUE)
colnames(mono@meta.data)[(length(colnames(mono@meta.data))-(length(scorelist)-1)):length(colnames(mono@meta.data))]=names(scorelist)
colnames(mono@meta.data)
meta=mono@meta.data
export(meta,paste0(filep,'mono-score_meta.xlsx'))

meta=import(paste0(filep,'mono-score_meta.xlsx'))
meta=meta[,c('celltype','Pathogen_recognition_and_clearance','Pro_inflammatory','Antigen_presentation','MHC_II')]

meta_avg <- meta %>%
  group_by(celltype) %>%
  summarise(across(
    c(
      Pathogen_recognition_and_clearance, 
      Pro_inflammatory, 
      Antigen_presentation, 
      MHC_II
    ),
    mean, na.rm = TRUE
  ))

meta_mat <- as.data.frame(meta_avg)
rownames(meta_mat) <- meta_mat$celltype
meta_mat$celltype <- NULL
dfm=t(meta_mat)
rownames(dfm)=gsub('_',' ',rownames(dfm))
p2=pheatmap(dfm,
            cluster_rows = FALSE,
            color = colorRampPalette(c('#8ccbea','#f7df87','#e68b81'))(100),
            cluster_cols = FALSE,
            scale = "row",  
            fontsize_row = 10,
            fontsize_col = 10,
            main = "")

## sFigure 4D
Idents(sce.all)=sce.all$celltype
nk1=subset(sce.all,idents=c("NK_c01-CD56high"))
nk2=subset(sce.all,idents=c("NK_c02-CD56low"))
nk1$celltype=as.vector(nk1$celltype)
nk2$celltype=as.vector(nk2$celltype)

genes=c("IFNG", 'TNF', "CCL4", "CCL4L2", "CCL5", "IL32", "IL16", "CCL3",  
        "KIR3DL1", "KIR3DL2", "KIR2DL1", "KIR2DL3", "KIR2DL1", "KIR2DL1")
average_expression <- AggregateExpression(nk1,
                                          features =genes,
                                          assays = 'RNA',
                                          group.by = c('outcome'),
                                          return.seurat = FALSE)
average_expression <- average_expression$RNA

expression_matrix <- as.matrix(average_expression)
hm1=expression_matrix[,c(2,1,3)]
hm1=hm1[genes,]
pheatmap(hm1,
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              color = colorRampPalette(c('#8ccbea','#f7df87','#e68b81'))(100),
              angle_col ='0',
              gaps_row = c(8),
              scale='row',
              main = 'NK_c01-CD56high')

genes=c("NKG7", "PRF1", "GSDMD", "GZMB", "GZMA", "FASLG", "TNFSF10", "KLRK1", 'FCGR3A', 'ITGAL', 'ITGB2')
average_expression <- AggregateExpression(nk2,
                                          features =genes,
                                          assays = 'RNA',
                                          group.by = c('outcome'),
                                          return.seurat = FALSE)
average_expression <- average_expression$RNA
expression_matrix <- as.matrix(average_expression)
hm2=expression_matrix[,c(2,1,3)]
hm2=hm2[genes,]
pheatmap(hm2,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               color = colorRampPalette(c('#8ccbea','#f7df87','#e68b81'))(100),
               angle_col ='0',
               #gaps_row = c(7,11),
               scale='row',
               main = 'NK_c02-CD56low')

## sFigure 4E
Idents(sce.all)=sce.all$celltype_major
dc=subset(sce.all,idents = c("DC"))
dc$outcome=factor(dc$outcome,levels = c('Control','Alive','Deceased'))

pdc=DotPlot(dc,features =unique( 
  c(
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", 
    "HLA-DQB1", 
    'CD83', 'IRF8'
  )),
  group.by = "outcome")
pdc

data1 <- pdc$data
data1$Average_Expression <- scales::rescale(data1$avg.exp.scaled, to = c(-1, 1))
data1$Percent_Expressed <- scales::rescale(data1$pct.exp, to = c(0.25, 1))

ggplot(data1, aes(x = features.plot, y = id, size = Percent_Expressed, color = Average_Expression)) +
  geom_point() +
  scale_size(name = "Percent_Expressed", range = c(2, 5)) +
  scale_color_gradient(name = "Average_Expression", low = '#b5dbfd', high ='#f7df87') +
  coord_flip() +
  labs(x = '', y = '', title = "mDC") +
  theme_classic() +
  th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))