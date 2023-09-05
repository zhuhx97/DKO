
rm(list = ls())
######*    install   *######
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(DESeq2)
library(dplyr)
source("/mnt/HDDN8R6/zhuhaoxian/RNAseq_raw/CD8/myvolcano.R")
source("/mnt/HDDN8R6/zhuhaoxian/RNAseq_raw/CD8/RNAseq_diff_analysis.R")

######***********************######
######*   1.data load  *######
######***********************######
data.matrix <- read.table("/mnt/HDDN8R6/zhuhaoxian/RNAseq_raw/LIVER/all_rpkm/merge_rpkm/merged_rpkm.txt"
                       ,sep = "\t")  
colnames(data.matrix)=data.matrix[1,] 
data.matrix=data.matrix[-c(1),]

rownames(data.matrix)=data.matrix[,1]
data.matrix=data.matrix[,-c(1)]

class(data.matrix)
for (i in 1:ncol(data.matrix)) {
  data.matrix[,i]=as.numeric(data.matrix[,i]) 
}

data.matrix.scale=as.data.frame(t(scale(t(data.matrix)))) 
data.matrix.scale=na.omit(data.matrix.scale)
data.matrix.scale.t=t(data.matrix.scale)



######**************************#####
######*   3.group analysis  *######
######*************************######

## counts data ----

data.count <- read.table("/mnt/HDDN8R6/zhuhaoxian/RNAseq_raw/LIVER/all_count/merge_count/merge_count.txt"
                          ,sep = "\t") 
colnames(data.count)=data.count[1,]
data.count=data.count[-c(1),]

rownames(data.count)=data.count[,1] 
data.count=data.count[,-c(1)] 


for (i in 1:ncol(data.count)) {
  data.count[,i]=as.numeric(data.count[,i]) 
} 
data.count=data.count[,-c(11:17)] 

## DKO vs Ctrl ----
dd.group = data.frame(Sample = factor(c("CTRL",rep("CART",5),rep("CTRL",4)), 
                                       levels = c("CTRL", "CART")),
                       row.names = colnames(data.count)) 
data.diff=RNAseq_diff(count.mtx = data.count,group.info = dd.group,
                      count.cutoff = 0,volcano.geneno=100,volcano.textsize=0.2) 
write.table(data.diff,file = "DEG_CTRL_vs_CART.txt",
            sep = "\t",col.names = NA)

data.diff.sig=subset(data.diff,subset=significant!="unchanged")
data.diff.sig=data.diff.sig[order(data.diff.sig$log2FoldChange,decreasing = T),]
up_diff_matrix=data.diff.sig[data.diff.sig$log2FoldChange>=1,] 
down_diff_matrix=data.diff.sig[data.diff.sig$log2FoldChange<(-1),] 

upgene.name=rownames(up_diff_matrix)
downgene.name=rownames(down_diff_matrix)
heatmap.genes=c(rownames(head(data.diff.sig,20)),
                rownames(tail(data.diff.sig,20))) 
heatmap.genes=c("Acta2","Col1a1","Timp1","Tgfb1")
data.matrix.scale.clean <- data.matrix.scale[,-c(11:17)]
                         
heatmap.matrix=data.matrix.scale.clean[heatmap.genes,] 

pdf("DEG_heatmap.pdf", width = 8, height = 8) 
pheatmap(heatmap.matrix,cluster_rows = T,cluster_cols = F,
         border_color = "black",treeheight_row = 0,
         color =rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))
dev.off()

######**************************#####
######*        4.GO分析       *######
######*************************######
require("org.Mm.eg.db")
require(DOSE)
require(clusterProfiler)
species="mouse" 
if(species=="h"|species=="human"){
  organism_db="org.Hs.eg.db"
  organism="hsa"
}else{
  organism_db="org.Mm.eg.db"
  organism="mmu"
}
up_geneid = bitr(upgene.name, fromType="SYMBOL", 
                 toType="ENTREZID", OrgDb=organism_db)
up_geneid = up_geneid$ENTREZID 
up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db, 
                   keyType = 'ENTREZID', ont="BP",
                   pvalueCutoff=0.05,readable=TRUE) 
down_geneid = bitr(downgene.name, fromType="SYMBOL", 
                   toType="ENTREZID", OrgDb=organism_db)
down_geneid = down_geneid$ENTREZID
down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db, 
                     keyType = 'ENTREZID', ont="BP",pvalueCutoff=0.05,readable=TRUE)
write.csv(up_ego,"GO_ConA_UP_result.csv") 
write.csv(down_ego,"GO_ConA_DOWN_result.csv")

pdf("GO_UP_top20_barplot.pdf", width =8, height = 8)
barplot(up_ego, showCategory=20,title=paste0("GO_UP"))
dev.off()
pdf("GO_UP_top20_dotplot.pdf", width =8, height = 8)
dotplot(up_ego,showCategory=10,title=paste0("GO_UP")) 
dev.off()
pdf("GO_down_top20_barplot.pdf", width =8, height = 8)
barplot(down_ego, showCategory=20,title=paste0("GO_DOWN"))
dev.off()
pdf("GO_down_top20_dotplot.pdf", width =8, height = 8)
dotplot(down_ego,showCategory=15,title=paste0("GO_DOWN"))+
  scale_colour_gradientn(colours =rev(color1(1000)))
dev.off()

######***********************######
######*    GSEA analysis    *######
######***********************######
write.table(data.matrix,"gsea_exprdata.txt",sep="\t",col.names = NA)
dim(data.matrix)

######***********************######
######*    GSEA analysis    *######
######***********************######
gsea_score <- read.table('/mnt/HDDN8R6/zhuhaoxian/RNAseq_raw/LIVER/result/gsea_score.txt',header = T,sep = '\t')
head(gsea_score)
gsea_score$NOM.p.val[gsea_score$NOM.p.val==0] <- 1e-4
color1 <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdBu")))

gsea_score$log10 <- -log10(gsea_score$NOM.p.val)
ggplot(gsea_score,aes(x = factor(gsea_score$celltype,levels = c('CTRL','CAR')), y = factor(gsea_score$NAME,levels = unique(gsea_score$NAME)[rev(c(2,6,4,5,1,3))])
                           , color = NES,size = log10, fill = NES)) + 
  geom_point() + theme_bw() + 
  # scale_color_(low = 'red',high = 'blue') +
  #scale_size_area(max_size=20) +
  scale_size_continuous(range = c(0.5,10)) + 
  xlab('celltype') + ylab('') + 
  #scale_colour_gradient(low = 'blue', high = 'red') + 
  scale_colour_gradientn(colours =color1(1000)) +  
  coord_equal(ratio = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 


