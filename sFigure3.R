library(CellChat)
filep="./MiddleFile/"

ctl_big=c('B','Plasma','CD4+T','CD8+T','NK','DC','Mono')
ctc_big=c('#f091a0',"#d35b7e",'#fab378','#80b5b8','#caadd8','#80b1d2','#6B4C3D')
names(ctc_big)=ctl_big
ctcol=ctc_big

pt=c("Control","Alive","Deceased")
hc=readRDS(paste0(filep,pt[1],"-PBMC-cellchat.rds"))
al=readRDS(paste0(filep,pt[2],"-PBMC-cellchat.rds"))
de=readRDS(paste0(filep,pt[3],"-PBMC-cellchat.rds"))

hc <- aggregateNet(hc)
al <- aggregateNet(al)
de <- aggregateNet(de)

##sFigure 3A
source("./new_netVisual_heatmap.R")
r=c(0, 0.02)
mr=0.08
mc=0.08

hc_heatmap <- new_netVisual_heatmap(hc, measure = "weight", title.name = "HC", 
                                    font.size = 10, font.size.title = 10, 
                                    max_row_barplot =mr ,
                                    max_col_barplot=mc,
                                    color.use = ctcol,
                                    color.heatmap = cht, color_range = r)
al_heatmap <- new_netVisual_heatmap(al, measure = "weight", title.name = "Alive", 
                                    font.size = 10, font.size.title = 10, 
                                    max_row_barplot =mr ,
                                    max_col_barplot=mc,
                                    color.use = ctcol,
                                    color.heatmap = cht, color_range = r)
de_heatmap <- new_netVisual_heatmap(de, measure = "weight", title.name = "Deceased", 
                                    font.size = 10, font.size.title = 10, 
                                    max_row_barplot =mr ,
                                    max_col_barplot=mc,
                                    color.use = ctcol,
                                    color.heatmap = cht, color_range = r)

##sFigure 3B
netVisual_aggregate(hc, signaling = 'APRIL', 
                    color.use =ctcol,
                    layout = "chord")

##sFigure 3C
netVisual_aggregate(de, signaling = 'APRIL', 
                    color.use =ctcol,
                    layout = "chord")

##sFigure 3D-E
p1=netVisual_bubble(hc, 
                    signaling ='APRIL',
                    sources.use = c(6,7), 
                    targets.use = c(2), remove.isolate = FALSE)
p2=netVisual_bubble(de, 
                    signaling ='APRIL',
                    sources.use = c(6,7), 
                    targets.use = c(2), remove.isolate = FALSE)
p=p1+p2
p
