percentplot_sample=function(seurat_object,filename,width=6,height=4,
                            legend.key.size=4,sample_name=NULL,flip=F,angle=45,
                            colors=colors2){
  num <- table(seurat_object@active.ident)
  if(is.null(sample_name)){
    sample_name=unique(seurat_object@meta.data$samples)
  }
  #sample_name=c("C1","C9","C10","Others")
  frequency_matrix=data.frame(names(num))
  cluster=levels(seurat_object)
  for (i in sample_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@meta.data$samples==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@active.ident)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(levels(seurat_object@active.ident),
                 names(table(subset_data@active.ident)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    row.names(frequency)=frequency$Var1
    frequency=frequency[cluster,]
    frequency_matrix[,i]=frequency$Freq
  }
  names(frequency_matrix)=c("cluster_name",sample_name)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  colnames(frequency_matrix_m)=c("Cluster","Sample","Frequency","count","Percent")
  #colorm=cbind(names(num),colors[names(num)])
  colorm=cbind(names(num),colors[1:length(names(num))])
  row.names(colorm)=colorm[,1]
  color_percent=as.character(colorm[,2])
  frequency_matrix_m2=as.data.frame(frequency_matrix_m)
  p <- ggplot(frequency_matrix_m2, aes(x=Sample, y=Percent, group=Cluster)) +
    geom_bar(stat="identity", color="white",position="fill", aes(fill=Cluster))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = angle),
          legend.key.size=unit(legend.key.size,'mm'))
  if(flip==TRUE){
    p=p+coord_flip()
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height*2)
  print(p)
  dev.off()
}

percentplot_cluster_group=function(seurat_object,filename,width=6,height=4,flip=F){
  num <- table(seurat_object@meta.data$group)
  cluster_name=levels(seurat_object)
  #cluster_name=sort(cluster_name)
  frequency_matrix=data.frame(names(num))
  group=frequency_matrix$names.num.
  for (i in cluster_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@active.ident==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@meta.data$group)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(as.character(unique(seurat_object$group)),
                 names(table(subset_data$group)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    suppressWarnings(if(!is.na(as.numeric(as.character(frequency$Var1)[1]))){
      print("reorder cluster name")
      frequency=frequency[order(as.numeric(as.character(frequency$Var1))),]
    })
    row.names(frequency)=frequency$Var1
    frequency=frequency[group,]
    frequency_matrix[,i]=frequency$Freq
  }
  cluster_name2=as.character(cluster_name)
  names(frequency_matrix)=c("cluster_name",cluster_name2)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  #frequency_matrix_m=frequency_matrix_m[sort(frequency_matrix_m$Clusters),]
  colnames(frequency_matrix_m)=c("Group","Cluster","Frequency","count","Percent")
  colorm=cbind(names(num),color_blue[5:length(names(num))])
  row.names(colorm)=colorm[,1]
  colorm=colorm[sort(names(num)),]
  color_percent=as.character(colorm[,2])
  p <- ggplot(frequency_matrix_m, aes(x=Cluster, y=Percent, group=Group)) +
    geom_bar(stat="identity", color="white",position="fill", aes(fill=Group))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 45))
  if(flip==TRUE){
    p=p+coord_flip()
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height)
  print(p)
  dev.off()
}

percentplot_cluster_sample=function(seurat_object,filename,width=6,height=4,flip=F){
  num <- table(seurat_object@meta.data$samples)
  cluster_name=levels(seurat_object)
  sample_name=unique(seurat_object$samples)
  #cluster_name=sort(cluster_name)
  #sample_name=c(non_mut,NOTCH_mut)
  frequency_matrix=data.frame(names(num))
  for (i in cluster_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@active.ident==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@meta.data$samples)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(as.character(unique(seurat_object$samples)),
                 names(table(subset_data$samples)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    suppressWarnings(if(!is.na(as.numeric(as.character(frequency$Var1)[1]))){
      print("reorder cluster name")
      frequency=frequency[order(as.numeric(as.character(frequency$Var1))),]
    })
    row.names(frequency)=frequency$Var1
    frequency=frequency[sample_name,]
    frequency_matrix[,i]=frequency$Freq
  }
  cluster_name2=as.character(cluster_name)
  names(frequency_matrix)=c("cluster_name",cluster_name2)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  colnames(frequency_matrix_m)=c("Sample","Cluster","Frequency","count","Percent")
  colorm=cbind(names(num),c(color_pink,color_yellow)[1:length(names(num))])
  row.names(colorm)=colorm[,1]
  #colorm=colorm[sort(names(num)),]
  color_percent=as.character(colorm[,2])
  frequency_matrix_m2=as.data.frame(frequency_matrix_m)
  p <- ggplot(frequency_matrix_m2, aes(x=Cluster, y=Percent, group=Sample)) +
    geom_bar(stat="identity", position="fill",color="white",aes(fill=Sample))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 45))
  if(flip==TRUE){
    levels(frequency_matrix_m2$Cluster)=rev(levels(frequency_matrix_m2$Cluster))
    p <- ggplot(frequency_matrix_m2, aes(x=rev(Cluster), y=Percent, group=Sample)) +
      geom_bar(stat="identity", position="fill",color="white",aes(fill=Sample))+
      scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
            axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 45))+
      coord_flip()+xlab("Cluster")
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height)
  print(p)
  dev.off()
}

percent_plot_wangfei_group=function(seurat_object,filename){
  for (j in 23:31){
    gn=names(seurat_object@meta.data)[j]
    print(gn)
    num <- table(seurat_object@active.ident)
    #sample_name=unique(seurat_object@meta.data$samples)
    group_name=unique(seurat_object@meta.data[,j])
    frequency_matrix=data.frame(names(num))
    for (i in group_name) {
      cell_name=row.names(seurat_object@meta.data[seurat_object@meta.data[,j]==i,])
      subset_data=subset(seurat_object, cells = cell_name)
      frequency=table(subset_data@active.ident)/ncol(subset_data)
      frequency
      frequency=as.data.frame(frequency)
      zero=setdiff(as.character(unique(seurat_object@active.ident)),
                   names(table(subset_data@active.ident)))
      if(length(zero)!=0){
        zero=as.data.frame(zero)
        zero$Freq=c(rep(0,nrow(zero)))
        names(zero)[1]="Var1"
        frequency=rbind(frequency,zero)
        print("reorder cluster name")
        frequency=frequency[order(as.numeric(as.character(frequency$Var1))),]
      }
      frequency=frequency[,-1]
      frequency_matrix=cbind(frequency_matrix,frequency)
    }
    names(frequency_matrix)=c("cluster_name",group_name)
    row.names(frequency_matrix)=frequency_matrix[,1]
    write.csv(frequency_matrix,paste0(filename,gn,"_frequency_matrix.csv"),row.names = F)
    library(reshape2)
    data_rownames <- as.character(names(num))
    data_colnames <- colnames(frequency_matrix)
    frequency_matrix$cluster <- data_rownames
    frequency_matrix=frequency_matrix[,-1]
    frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
    frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
    frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
      dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
    colnames(frequency_matrix_m)=c("Cluster","Group","Frequency","count","Percent")
    colorm=cbind(names(num),colors[1:length(names(num))])
    row.names(colorm)=colorm[,1]
    colorm=colorm[sort(names(num)),]
    color_percent=as.character(colorm[,2])
    p <- ggplot(frequency_matrix_m, aes(x=Group, y=Percent, group=Cluster)) +
      geom_bar(stat="identity", position="fill", aes(fill=Cluster))+
      scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
            axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 45))
    four=c("stage","differentiation")
    if(gn%in%four){
      pdf(paste0("percent_plot_",filename,gn,".pdf"),width =6, height = 4)
      print(p)
      dev.off()
    }else{
      pdf(paste0("percent_plot_",filename,gn,".pdf"),width =3, height = 4)
      print(p)
      dev.off()
    }
  }
}

percentplot_group=function(seurat_object,filename,width=6,height=4,
                            legend.key.size=4,sample_name=NULL,flip=F,angle=45,
                            colors=colors2){
  num <- table(seurat_object@active.ident)
  if(is.null(sample_name)){
    sample_name=unique(seurat_object@meta.data$group)
  }
  #sample_name=c("Normal","PE")
  frequency_matrix=data.frame(names(num))
  cluster=levels(seurat_object)
  for (i in sample_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@meta.data$group==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@active.ident)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(levels(seurat_object@active.ident),
                 names(table(subset_data@active.ident)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    row.names(frequency)=frequency$Var1
    frequency=frequency[cluster,]
    frequency_matrix[,i]=frequency$Freq
  }
  names(frequency_matrix)=c("cluster_name",sample_name)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  colnames(frequency_matrix_m)=c("Cluster","Sample","Frequency","count","Percent")
  #colorm=cbind(names(num),colors[names(num)])
  colorm=cbind(names(num),colors[1:length(names(num))])
  row.names(colorm)=colorm[,1]
  color_percent=as.character(colorm[,2])
  frequency_matrix_m2=as.data.frame(frequency_matrix_m)
  p <- ggplot(frequency_matrix_m2, aes(x=Sample, y=Percent, group=Cluster)) +
    geom_bar(stat="identity", color="white",position="fill", aes(fill=Cluster))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = angle),
          legend.key.size=unit(legend.key.size,'mm'))
  if(flip==TRUE){
    p=p+coord_flip()
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height*2)
  print(p)
  dev.off()
}


percentplot_choose=function(seurat_object,filename,width=6,height=4,
                           legend.key.size=4,sample_name=NULL,flip=F,angle=45,
                           colors=colors2,choose=NULL){
  num <- table(seurat_object@active.ident)
  if(is.null(sample_name)){
    sample_name=levels(tnk.combined@meta.data[,choose])
  }
  #sample_name=c("Normal","PE")
  frequency_matrix=data.frame(names(num))
  cluster=levels(seurat_object)
  for (i in sample_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@meta.data[,choose]==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@active.ident)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(levels(seurat_object@active.ident),
                 names(table(subset_data@active.ident)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    row.names(frequency)=frequency$Var1
    frequency=frequency[cluster,]
    frequency_matrix[,i]=frequency$Freq
  }
  names(frequency_matrix)=c("cluster_name",sample_name)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  colnames(frequency_matrix_m)=c("Cluster","Sample","Frequency","count","Percent")
  #colorm=cbind(names(num),colors[names(num)])
  colorm=cbind(names(num),colors[1:length(names(num))])
  row.names(colorm)=colorm[,1]
  color_percent=as.character(colorm[,2])
  frequency_matrix_m2=as.data.frame(frequency_matrix_m)
  p <- ggplot(frequency_matrix_m2, aes(x=Sample, y=Percent, group=Cluster)) +
    geom_bar(stat="identity", color="white",position="fill", aes(fill=Cluster))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = angle),
          legend.key.size=unit(legend.key.size,'mm'))
  if(flip==TRUE){
    p=p+coord_flip()
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height*2)
  print(p)
  dev.off()
}








