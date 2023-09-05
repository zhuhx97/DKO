setwd('C:/Users/dell/Desktop/DKO SC')
#引入所需R包
install.packages("BiocManager")
BiocManager::install("g++")
BiocManager::install('monocle')
BiocManager::install('export')
BiocManager::install("Seurat")
BiocManager::install("GSEABase")
BiocManager::install("zoo")
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(homologene)
library(GSEABase)
source('./Get_Monocle2.R')
###################

minFeatureRNA <- 200
maxFeatureRNA <- 4000
Percent.mt <- 0.5

E0_1_object <- subset(E0_1_object, subset = ((nFeature_RNA > minFeatureRNA) & (nFeature_RNA < maxFeatureRNA) & (percent.mt < Percent.mt)) )
print(dim(E0_1_object))

E0_2_object[["percent.mt"]] <- PercentageFeatureSet(E0_2_object, pattern = "^Mt")
VlnPlot(E0_2_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E0_2_object <- subset(E0_2_object, subset = ((nFeature_RNA > minFeatureRNA) & (nFeature_RNA < maxFeatureRNA) & (percent.mt < Percent.mt)) )

library(dplyr)
#E0_1
E0_1_object=E0_1_object%>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%RunPCA()
ElbowPlot(object = E0_1_object,ndims = 50)
E0_1_object <- FindNeighbors(object = E0_1_object, dims = 1:30)
E0_1_object <- FindClusters(object = E0_1_object, resolution = 0.4)

E0_1_object <- RunUMAP(E0_1_object, reduction = "pca", dims = 1:30)
E0_1_object <- RunTSNE(object = E0_1_object, dims = 1:30)
DimPlot(E0_1_object,reduction = 'umap',label = T)

#E0_2
E0_2_object=E0_2_object%>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%RunPCA()
ElbowPlot(object = E0_2_object,ndims = 50)
E0_2_object <- FindNeighbors(object = E0_2_object, dims = 1:30)
E0_2_object <- FindClusters(object = E0_2_object, resolution = 0.4)
#run tsne and umap
E0_2_object <- RunUMAP(E0_2_object, reduction = "pca", dims = 1:30)
E0_2_object <- RunTSNE(object = E0_2_object, dims = 1:30)
DimPlot(E0_2_object,reduction = 'umap',label = T)
source("/mnt/HDDN8R6/zhuhaoxian/R/Get_integrated.R")
sce.big=Get_integrated(list(E0_1_object,E0_2_object))

CD8T <- subset(sce.big,idents = c(0,10,7,1,2,3,8,9))
CD8T <- ScaleData(CD8T,features = rownames(CD8T))
CD8T <- ScaleData(CD8T,features = rownames(CD8T))
CD8T <- RunPCA(CD8T) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.15)
CD8T <- FindClusters(CD8T, resolution = 0.6)
CD8T<- RunUMAP(CD8T, dims = 1:30)
CD8T <- RunTSNE(CD8T, dims = 1:30)
myfindmarkers(DKOCD8T_clean,filename="AllCD8subsetnew",colors = colors)
table(CD8T$samples)
DimPlot(CD8T,reduction = 'tsne',split.by = 'samples',cols = colors,pt.size = 1,label = T)
DimPlot(CD8T,reduction = 'umap',cols = colors,label = F,pt.size = 1)
FeaturePlot(CD8T,reduction = 'tsne',features = c('Cd8a','Cd8b1','Mki67','Cd44','Sell','S1pr1','Ccr7','Ly6c2','Cxcr6','Runx3','Klf2','Tox'),
            cols = c('lightgrey','#F4A582',"#D6604D",'#B2182B'),ncol = 4)
FeaturePlot(CD8T,reduction = 'tsne',features = c('Cd8a','Cd8b1','Top2a','Mki67'),cols = c('lightgrey','#F4A582',"#D6604D",'#B2182B'))
FeaturePlot(CD8T,reduction = 'tsne',features = c( 'Pdcd1','Gzmb','Gzma','Ifng','Prf1','Fasl'))
percentplot_sample(CD8T,filename="percent_CD8_WTvsDKO",colors = colors,width = 4,height = 2)



FeaturePlot(CD8T,reduction = 'tsne',features = c("Cd8b1",'Gzma','Prf1','Ccl5','Rgs1','Cxcr6','Itgal','Ly6c2','Il7r','Sell','Cx3cr1','Mki67','Top2a','Tcf7','Ccr7'),
            cols = c('lightgrey','#F4A582',"#D6604D",'#B2182B'),ncol = 5)
source("C:/Users/dell/Desktop/myfindmarkers.R")
myfindmarkers(CD8T,filename="AllCD8subsetnew",colors = colors)

library(monocle)
DKOcd8t.cds <- Get_monocle2(CD8T)
DKOcd8t.cds <- orderCells(DKOcd8t.cds,reverse = T)
plot_cell_trajectory(DKOcd8t.cds, color_by = "Pseudotime",cell_size = 2)

plot_cell_trajectory(DKOcd8t.cds, color_by = "RNA_snn_res.0.6",cell_size = 2) + 
  scale_color_manual(values = colors2) 

plot_cell_trajectory(DKOcd8t.cds, color_by = "RNA_snn_res.0.6",cell_size = 2) + 
  scale_color_manual(values = colors2)  +  facet_wrap('RNA_snn_res.0.6',ncol = 3)
brewer.pal(8,'Set2')
dimplot(DKO_CD8,)
#单独DKO的Trm分析
DKO_CD8Trm <- subset(CD8T,idents='0')
DimPlot(DKO_CD8Trm,reduction = "umap",pt.size = 1,cols = colors2[c(2,7,8)])
DKO_CD8Trm <- NormalizeData(DKO_CD8Trm, normalization.method = "LogNormalize", scale.factor = 10000)

DKO_CD8Trm <- FindVariableFeatures(DKO_CD8Trm, selection.method = "vst", nfeatures = 2000)
# data.seurat <- CellCycleScoring(data.seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DKO_CD8Trm <- ScaleData(DKO_CD8Trm,features = rownames(DKO_CD8Trm))
DKO_CD8Trm <- RunPCA(DKO_CD8Trm) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.6)
DKO_CD8Trm <- FindClusters(DKO_CD8Trm, resolution = 0.3)
DKO_CD8Trm <- RunUMAP(DKO_CD8Trm, dims = 1:30)
DKO_CD8Trm <- RunTSNE(DKO_CD8Trm, dims = 1:30)
DimPlot(DKO_CD8Trm,cols =colors,pt.size = 2)
FeaturePlot(DKO_CD8Trm,reduction = 'umap', features = c('Cd200r1','Cd244a'),pt.size = 1)
VlnPlot(DKO_CD8Trm,features =c('Cd244a','Cd200r1'),
        pt.size =0,
        cols = colors2[c(2,7,8)],ncol=1)


for (i in dir('./TO_YSH/gmt/trm/')) {
  
  gspath=paste0('./TO_YSH/gmt/trm/',i)
  geneset=getGmt(gspath)
  geneset[[1]]@geneIds
  
  CD8T <- AddModuleScore(CD8T,features = list(human2mouse(geneset[[1]]@geneIds)$mouseGene), name = geneset[[1]]@setName)
  
}
for (i in dir('./TO_YSH/gmt/trmdown/')) {
  
  gspath=paste0('./TO_YSH/gmt/trmdown/',i)
  geneset=getGmt(gspath)
  geneset[[1]]@geneIds
  
  CD8T <- AddModuleScore(CD8T,features = list(human2mouse(geneset[[1]]@geneIds)$mouseGene), name = geneset[[1]]@setName)
  
}

colnames(CD8T@meta.data)

FeaturePlot(CD8T,reduction = 'tsne',features = c('TRM_SIGNATURE_UP1'),pt.size = 1.5,
            cols = c('lightgrey','#F4A582',"#D6604D",'#B2182B') )
FeaturePlot(CD8T,reduction = 'tsne',features = c('TRM_SIGNATURE_DOWN1'),pt.size = 1.5,
            cols = c('lightgrey','#F4A582',"#D6604D",'#B2182B') )
##############
### heatmap ----
DKO_Trm.ave=AverageExpression(DKO_CD8Trm,return.seurat = T)
DKO_Trm.avetr=AverageExpression(DKO_CD8Trm,return.seurat = F)
DKO_Trm.ave$RNA[1:5,]
class(DKO_Trm.ave)
DKO_Trm.ave2tr=as.data.frame(DKO_Trm.avetr$RNA)
exhaust=c('Cd160','Lag3','Pdcd1','Cd244a','Tigit','Tox','Id2','Cd101')
stem=c('Tcf7','Slamf6','Id3','Eomes','Il7r','Cxcr5','Cxcr3','Bcl6')
cytotoxic=c('Gzmb','Gzmk','Fasl','Prf1','Ifng','Il2','Ccl5','Nkg7')
costimu=c('Icos','Cd28','Tnfsf14','Tnfrsf4')

DKO_Trm.ave2trexhaust=DKO_Trm.ave2tr[exhaust,]
DKO_Trm.ave2trstem=DKO_Trm.ave2tr[stem,]
DKO_Trm.ave2trcytotoxic=DKO_Trm.ave2tr[cytotoxic,]
DKO_Trm.ave2trcostimu=DKO_Trm.ave2tr[costimu,]

DKO_Trm.ave2trexhaust=as.data.frame(t(scale(t(DKO_Trm.ave2trexhaust))))
DKO_Trm.ave2trstem=as.data.frame(t(scale(t(DKO_Trm.ave2trstem))))
DKO_Trm.ave2trcytotoxic=as.data.frame(t(scale(t(DKO_Trm.ave2trcytotoxic))))
DKO_Trm.ave2trcostimu=as.data.frame(t(scale(t(DKO_Trm.ave2trcostimu))))

pheatmap(DKO_Trm.ave2trexhaust,border_color = "black",
         treeheight_col = 0,treeheight_row = 0,cluster_cols = F,cellwidth = 20,cellheight = 20)
pheatmap(DKO_Trm.ave2trstem,border_color = "black",
         treeheight_col = 0,treeheight_row = 0,cluster_cols = F,cellwidth = 20,cellheight = 20)
pheatmap(DKO_Trm.ave2trcytotoxic,border_color = "black",
         treeheight_col = 0,treeheight_row = 0,cluster_cols = F,cellwidth = 20,cellheight = 20)
pheatmap(DKO_Trm.ave2trcostimu,border_color = "black",
         treeheight_col = 0,treeheight_row = 0,cluster_cols = F,cellwidth = 20,cellheight = 20)


DimPlot(CD8T,reduction = 'tsne',cols = colors,pt.size = 1,label = T)

source("./myvolcano.R")
CD8Tnew_volcano <- RenameIdents(CD8T,
                                '0'='Trm',
                                '1'='nonTrm',
                                '2'='nonTrm',
                                '3'='nonTrm',
                                '4'='nonTrm',
                                '5'='nonTrm')
DimPlot(CD8Tnew_volcano,reduction = 'tsne',cols = colors,pt.size = 1,label = T)

cell_20_gene=FindMarkers(CD8Tnew_volcano,ident.1 = "Trm",ident.2 = "nonTrm",logfc.threshold = 0.25)
cell_20_gene=cell_20_gene[order(cell_20_gene$avg_log2FC,decreasing = F),]
DotPlot (CD8Tnew_volcano, features= c("Cxcr6","Pdcd1","Lag3"),dot.scale = 8)+ scale_color_gradientn(values= seq(0,1,0.2), colours = c('#330066','#336699','#66CC66','#FFCC33'))
FeaturePlot(CD8Tnew_volcano,reduction = 'tsne', features = c("Cxcr6","Pdcd1","Lag3"),pt.size = 1)
VlnPlot(CD8Tnew_volcano,features = c("Cxcr6","Pdcd1","Lag3"))
myvolcano(cell_20_gene,filename = '7-0-1',gene.plot = 0,p.cutoff=0.05,text.size = 0)
write.csv(cell_20_gene,file = '7-0.csv')
