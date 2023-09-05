#' To integrate different batch seurat.obj together
#'
#' @param umi.list the seurat list
#' @importFrom Seurat FindIntegrationAnchors IntegrateData
#' @return the integrated seurat.obj
#' @export
#'
#' @examples
#' \dontrun{
#' seurat.int = Get_integrated(list(seurat.obj1,seurat.obj2,seurat.obj3))
#' }
#'
Get_integrated = function(umi.list) {
  
  if (requireNamespace("future", quietly = TRUE)) {
    
    future::plan("multisession", workers = 1)
    umi <- FindIntegrationAnchors(object.list = umi.list, dims = 1:40)
    
    future::plan("multisession", workers = 15)
    umi.integrated <- IntegrateData(anchorset = umi, dims = 1:30)
    SeuratObject::DefaultAssay(umi.integrated) <- "integrated"
    future::plan("multisession", workers = 1)
    umi.integrated <- ScaleData(umi.integrated, verbose = T)
    future::plan("multisession", workers = 10)
    umi.integrated <- RunPCA(umi.integrated, verbose = T)
    
    DimPlot(object = umi.integrated, reduction = "pca")
    
    ElbowPlot(object = umi.integrated, ndims = 50)
    umi.integrated <- FindNeighbors(object = umi.integrated, dims = 1:30)
    umi.integrated <- FindClusters(object = umi.integrated, resolution = 0.6)
    colnames(umi.integrated@meta.data)
    
    ###tsne
    umi.integrated <- RunTSNE(object = umi.integrated, dims = 1:30,check_duplicates = FALSE)
    umi.integrated <- RunUMAP(umi.integrated, reduction = "pca", dims = 1:30)
    future::plan("multisession", workers = 1)
    
  }else{
    umi <- FindIntegrationAnchors(object.list = umi.list, dims = 1:40)
    
    umi.integrated <- IntegrateData(anchorset = umi, dims = 1:30)
    SeuratObject::DefaultAssay(umi.integrated) <- "integrated"
    umi.integrated <- ScaleData(umi.integrated, verbose = T)
    umi.integrated <- RunPCA(umi.integrated, verbose = T)
    
    DimPlot(object = umi.integrated, reduction = "pca")
    
    ElbowPlot(object = umi.integrated,ndims = 50)
    umi.integrated <- FindNeighbors(object = umi.integrated, dims = 1:30)
    umi.integrated <- FindClusters(object = umi.integrated, resolution = 0.6)
    
    ###tsne
    umi.integrated <- RunTSNE(object = umi.integrated, dims = 1:30,check_duplicates = FALSE)
    umi.integrated <- RunUMAP(umi.integrated, reduction = "pca", dims = 1:30)
  }
  
  return(umi.integrated)
  
}
