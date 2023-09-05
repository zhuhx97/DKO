#monocle 2-----
#BiocManager::install('monocle')
library(monocle)
seurat_object = CD8
num_expressed=1
dispersion_empiricals=0.1
mean_expressions=0.01
num_cells_expressed <- mean_expression <- dispersion_empirical <- NULL
###Monocle2 三件套
#表达矩阵
ct <- seurat_object@assays$RNA@data
sample_ann <- seurat_object@meta.data
gene_ann <- data.frame(
   "gene_short_name" = row.names(ct),
   row.names = row.names(ct))

###转换成AnnotationDataframe对象
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)

# 构建CDS对象
#expressionFamily ，选择数据分布:
#FPKM/TPM 值是log-正态分布的；UMIs和原始count值用负二项分布模拟的效果更好
#负二项分布有两种方法，这里选用了negbinomial.size，另外一种negbinomial稍微更准确一点，但速度大打折扣，它主要针对非常小的数据集
sc_cds <- newCellDataSet(
  as.matrix(ct),
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)

#####质控过滤
cds <- sc_cds

# 设置一个基因表达量的过滤阈值，结果会在cds@featureData@data中新增一列num_cells_expressed，记录这个基因在多少细胞中有表达
cds <- detectGenes(cds, min_expr = 0.1)
# 结果保存在cds@featureData@data
print(head(cds@featureData@data))
##基因过滤
#选出表达大于5的
expressed_genes <- row.names(subset(cds@featureData@data,
                                    num_cells_expressed >= num_expressed))
cds <- cds[expressed_genes,]
#####计算Size factors和"dispersion" values，
# Size factors可以帮助我们对不同细胞中mRNA表达的差异进行归一化处理
# "dispersion" values将有助于后续的基因差异表达分析。
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table <- dispersionTable(cds) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= mean_expressions) # 挑表达量不太低的
unsup_clustering_genes <- subset(unsup_clustering_genes, dispersion_empirical >= dispersion_empiricals) # 挑表达量不太低的

cluster_DEG_genes <- differentialGeneTest(cds[as.character(unsup_clustering_genes$gene_id),],fullModelFormulaStr = '~rename',cores = 6)
dim(cluster_DEG_genes)
cluster_DEG_gene <- rownames(cluster_DEG_genes)[order(cluster_DEG_genes$qval)][1:2000]

cds <- setOrderingFilter(cds, as.character(cluster_DEG_gene))  # 准备聚类基因名单
print(length(cluster_DEG_genes))
print(plot_ordering_genes(cds))

#if (!is.null(ordering_genes)) {
#  cds <- setOrderingFilter(cds, ordering_genes)
#}
###降维
# 默认使用DDRTree的方法
# 进行降维
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')

###细胞排序
cds <- orderCells(cds,reverse = F)
##可视化
pdf("monocle_plot_sample.pdf",width = 6.5,height = 5)
plot_cell_trajectory(cds, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(values = c("tomato2","darkslategray","royalblue","orange","grey"),breaks =  as.character(levels(seurat_object@active.ident))) +
  facet_wrap('seurat_clusters')

dev.off()
pdf("monocle_plot_total.pdf",width = 6.5,height = 5)
p <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  scale_color_manual(values = c("tomato2","darkslategray","royalblue","orange","grey"),breaks =  as.character(levels(seurat_object@active.ident))) #+
#facet_wrap('seurat_clusters')
print(p)
dev.off()

plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1.5,show_branch_points = FALSE)+
  scale_color_gradientn(values = seq(0,1,0.2),colours =  c('yellow2','orange','purple3','#330066'))
