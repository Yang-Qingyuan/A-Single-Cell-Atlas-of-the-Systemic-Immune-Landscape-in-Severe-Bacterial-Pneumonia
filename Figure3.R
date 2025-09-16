library(CellChat)
library(ggplot2)
library(NMF)
library(patchwork)

filep="./MiddleFile/"
th=theme(text = element_text(#family = "Arial",
  color="black"), 
  axis.title = element_text(size = 10,face = "bold"),
  plot.title = element_text(size=12,face = "bold"),
  plot.caption = element_text(size=8),
  axis.line = element_line(linewidth = 1.5),  
  axis.ticks = element_line(linewidth = 1) ,
  axis.text = element_text(size = 10, face = "bold", color = "black"))
outcomeCol=c('#7ac7e2','#fabb6e','#b7282e')
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

listt <- list(Control = hc, Alive = al,Deceased=de)
cellchatt <- mergeCellChat(listt, add.names = names(listt))

list1 <- list(Control = hc, Alive = al)
cellchat1 <- mergeCellChat(list1, add.names = names(list1))

list2 <- list(Alive = al,Deceased=de)
cellchat2 <- mergeCellChat(list2, add.names = names(list2))

## Figure 3A
p1 <- compareInteractions(cellchatt, show.legend = F, group = c(1,2,3), 
                          measure = "weight",
                          color.use = outcomeCol)+th
p1

## Figure 3B
netVisual_diffInteraction(cellchat1, weight.scale = T,
                          color.use = ctcol,
                          color.edge=rev(outcomeCol[1:2]),
                          measure = "weight")

## Figure 3D
netVisual_diffInteraction(cellchat2, weight.scale = T,
                          color.use = ctcol,
                          color.edge=rev(outcomeCol[2:3]), 
                          measure = "weight")

## Figure 3C, 3E
object.list=listt

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}

p3c=netAnalysis_diff_signalingRole_scatter(
  object.list,comparison = c(1,2),
  title ='Control vs. Alive' ,
  color.use =ctcol)+th
p3e=netAnalysis_diff_signalingRole_scatter(object.list,comparison = c(2,3),
                                           title ='Alive vs. Deceased' , 
                                           color.use =ctcol)+th
p3c
p3e

## Figure 3F (incoming)
  # HC
ptcol=c('#C5DFF4','#C9DCC4','#DAA87C','#F4EEAC')
names(ptcol)=paste0('Pattern ',1:4)

HC=hc
selectK(HC, pattern = "incoming")
# Select the pattern number where both the Cophenetic and Silhouette values start to suddenly decrease.
nPatterns = 4 
HC <- identifyCommunicationPatterns(HC, pattern = "incoming", k = nPatterns, 
                                    width = 5, height = 9, font.size = 6)


netAnalysis_river(HC, pattern = "incoming",
                  color.use = ctcol,
                  color.use.pattern =ptcol,
                  font.size = 4,
                  cutoff = 0.3,
                  font.size.title = 10)

  # Alive
selectK(al, pattern = "incoming")
# Select the pattern number where both the Cophenetic and Silhouette values start to suddenly decrease.
nPatterns = 3 
al <- identifyCommunicationPatterns(al, pattern = "incoming", k = nPatterns, 
                                    width = 5, height = 9, font.size = 6)
netAnalysis_river(al, pattern = "incoming",
                  color.use = ctcol,
                  color.use.pattern =ptcol,
                  font.size = 4,
                  cutoff = 0.38,
                  font.size.title = 10)

  # Deceased
selectK(de, pattern = "incoming")
# Select the pattern number where both the Cophenetic and Silhouette values start to suddenly decrease.
nPatterns = 2 
de <- identifyCommunicationPatterns(de, pattern = "incoming", k = nPatterns, 
                                    width = 5, height = 9, font.size = 6)
netAnalysis_river(de, pattern = "incoming",
                  color.use = ctcol,
                  color.use.pattern =ptcol,
                  font.size = 4,
                  font.size.title = 10)

## Figure 3G (outgoing)
  # HC
ptcol=c('#C5DFF4','#C9DCC4','#DAA87C','#F4EEAC')
names(ptcol)=paste0('Pattern ',1:4)

HC=hc
selectK(HC, pattern = "outgoing")
# Select the pattern number where both the Cophenetic and Silhouette values start to suddenly decrease.
nPatterns = 2 
HC <- identifyCommunicationPatterns(HC, pattern = "outgoing", k = nPatterns, 
                                    width = 5, height = 9, font.size = 6)
netAnalysis_river(HC, pattern = "outgoing",
                  color.use = ctcol,
                  color.use.pattern =ptcol,
                  font.size = 4,
                  font.size.title = 10)

  # Alive
selectK(al, pattern = "outgoing")
# Select the pattern number where both the Cophenetic and Silhouette values start to suddenly decrease.
nPatterns = 2 
al <- identifyCommunicationPatterns(al, pattern = "outgoing", k = nPatterns, 
                                    width = 5, height = 9, font.size = 6)
netAnalysis_river(al, pattern = "outgoing",
                  color.use = ctcol,
                  color.use.pattern =ptcol,
                  font.size = 4,
                  font.size.title = 10)

  # Deceased
selectK(de, pattern = "outgoing")
# Select the pattern number where both the Cophenetic and Silhouette values start to suddenly decrease.
nPatterns = 2
de <- identifyCommunicationPatterns(de, pattern = "outgoing", k = nPatterns, 
                                    width = 5, height = 9, font.size = 6)
netAnalysis_river(de, pattern = "outgoing",
                  color.use = ctcol,
                  color.use.pattern =ptcol,
                  font.size = 4,
                  font.size.title = 10)

## Figure 3H
gg1 <- rankNet(cellchat1, mode = "comparison", 
               measure = 'weight',
               color.use = outcomeCol[1:2],
               stacked = T, do.stat = TRUE)+th
gg2 <- rankNet(cellchat2, mode = "comparison",
               color.use = outcomeCol[2:3],
               measure = 'weight',
               stacked = T, do.stat = TRUE)+th
p3h=gg1+gg2
p3h

## Figure 3I
netVisual_aggregate(hc, signaling =  'IL16', 
                    color.use =ctcol,
                    layout = "chord")

netVisual_aggregate(al, signaling =  'IL16', 
                    color.use =ctcol,
                    layout = "chord")

netVisual_aggregate(de, signaling =  'IL16', 
                    color.use =ctcol,
                    layout = "chord")

## Figure 3J
netVisual_aggregate(hc, signaling =  'BAG', 
                    color.use =ctcol,
                    layout = "chord")

netVisual_aggregate(al, signaling =  'BAG', 
                    color.use =ctcol,
                    layout = "chord")

netVisual_aggregate(de, signaling =  'BAG', 
                    color.use =ctcol,
                    layout = "chord")

## Figure 3K-M
netVisual_aggregate(de, signaling = 'RESISTIN', 
                    color.use =ctcol,
                    layout = "chord")

netVisual_aggregate(de, signaling = 'IL1', 
                    color.use =ctcol,
                    layout = "chord")

netVisual_aggregate(de, signaling = 'GRN', 
                    color.use =ctcol,
                    layout = "chord")