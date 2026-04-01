# .libPaths(c('/share/apps/R/R-4.1.2-4/lib64/R/library_old'))  # 已注释，使用系统默认R库路径
suppressMessages({
library(dplyr)
library(plyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(SingleR)
library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="the seurat rds file")
argv <- add_argument(argv,"--cluster", help="the way of cluster,such as tsne or umap")
argv <- add_argument(argv,"--speices", help="the speices of project,such as Mouse or Human")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)

each_rds <- argv$rds
cluster <- argv$cluster
speices <- argv$speices
outdir <- argv$outdir

#color_panel<-c('#FFC1C1','#FFA500','#FFF68F','#7FFFD4','#B0E2FF','#FF83FA','#F08080','#54FF9F','#C0FF3E','#6495ED','#8470FF','#FF69B4',
#               '#FF1493','#FFDEAD','#0000FF','#8B1C62','#FF0000','#FFFF00','#00EE00','#B22222','#00BFFF','#006400','#000080','#CD5C5C',
#               '#9400D3','#836FFF','#008B8B','#CD5B45','#2E8B57','#473C8B','#BCD2EE','#363636','#F5F5DC','#53868B','#D8BFD8','#8B5A00',
#               '#8B658B','#C1CDC1','#EED5D2','#000000','#828282','#BC8F8F','#FF4040','#98F5FF','#90EE90')
color_panel <-c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
             '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
             '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
             '#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')
color_panel <-c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
             '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
             '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
             '#CCC9E6','#3A6963','#68A180','#61B934','#476D87','#9FA3A8','#968175','#BC8F8F','#D8BFD8')
#注释
each_rds<-readRDS(file = each_rds)

counts<-each_rds@assays$RNA@counts
ann<-each_rds@meta.data$orig.ident  ###Annotation of cell
clusters<-each_rds@meta.data$seurat_clusters ###
###
counts@Dimnames[[1]]<-gsub('^hg19-|^GRCh38-|^mm10-','',ignore.case = T,rownames(counts)) ###
#counts@Dimnames[[1]]<-gsub('^mm10-','',ignore.case = T,rownames(counts))
counts@Dimnames[[1]]<-gsub('^MT-|[-_\\s]{1,3}mt','',ignore.case = T,rownames(counts))
clusters<-each_rds@meta.data$seurat_clusters ###
Cellcluster_Anno = CreateSinglerObject(counts, annot = ann, 'YR', min.genes = 0,
                                       species = speices, citation = "",
                                       ref.list = list(), normalize.gene.length = F, variable.genes = 'de',
                                       fine.tune = F, do.signatures = F, clusters = clusters, do.main.types = T,
                                       reduce.file.size = T, numCores = 10)
###singler1
clusterAnnotaion=Cellcluster_Anno$singler[[1]]$SingleR.clusters.main$labels #
clusterAnn_df<-data.frame(Annotation=clusterAnnotaion,stringsAsFactors = F)
clusterAnn_df$ID<-rownames(clusterAnn_df) ###
row_length<-nrow(clusterAnnotaion)-1
rownames(clusterAnnotaion)<-paste('cluster',0:row_length,sep='')
clusterAnnotaion<-data.frame(ClusterID=rownames(clusterAnnotaion),Annotation=clusterAnnotaion[,1])
write.table(clusterAnnotaion,file=paste0(outdir,"/Cluster_Cells_Annotation_HPCA.txt"),quote=F,sep="\t",col.names=T,row.names=F) #save annotation of cluster
singlecellanno<-data.frame(Cell_ID=rownames(Cellcluster_Anno$singler[[1]]$SingleR.single.main$labels),Annotation=Cellcluster_Anno$singler[[1]]$SingleR.single.main$labels[,1], stringsAsFactors = F)##extract cell annotation
singlecellanno$cluster<-mapvalues(singlecellanno$Cell_ID,rownames(each_rds@meta.data),as.character(each_rds@meta.data$seurat_clusters))
write.table(singlecellanno,file=paste0(outdir,"/SinglecellAnnotation_HPCA.txt"),quote=F,sep="\t",col.names=F,row.names=F) #save cell annotation
###match
each_rds@meta.data$celltype1<-clusterAnn_df$Annotation[match(each_rds@meta.data$seurat_clusters,clusterAnn_df$ID)]
each_rds@meta.data$cluster_celltype1<-paste(each_rds@meta.data$seurat_clusters,each_rds@meta.data$celltype1,sep=':')
each_rds@meta.data$cluster_celltype1<-factor(each_rds@meta.data$cluster_celltype1,levels=unique(each_rds@meta.data[order(each_rds@meta.data$seurat_clusters),]$cluster_celltype1))
###cluster umap/tsne
###single cell umap/tsne
if(cluster=='tsne'){
    add_lable_plot1<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('celltype1'),label = T,raster = F)
    add_lable_plot2<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('cluster_celltype1'),label = T,raster = F)
    ggsave(filename=paste0(outdir,'/Annotation_tsne_HPCA.png'),add_lable_plot1,width=12,height=10)
    ggsave(filename=paste0(outdir,'/Annotation_cluster_tsne_HPCA.png'),add_lable_plot2,width=12,height=10)
    out_plot_single = SingleR.PlotTsne(Cellcluster_Anno$singler[[1]]$SingleR.single,each_rds@reductions$tsne@cell.embeddings, do.labels = F,do.letters = T)
    ggsave(filename=paste0(outdir,'/Single_cell_annot_tsne_HPCA.png'),out_plot_single$p,dpi=400,width=10,height=8)
    #ggsave(filename=paste0(outdir,'/Single_cell_annot_tsne_HPCA.pdf'),out_plot_single$p,dpi=400,width=10,height=8)
}
if(cluster=='umap'){
    add_lable_plot1<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('celltype1'),label = T,raster = F)
    add_lable_plot2<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('cluster_celltype1'),label = T,raster = F)
    ggsave(filename=paste0(outdir,'/Annotation_umap_HPCA.png'),add_lable_plot1,width=12,height=10)
    ggsave(filename=paste0(outdir,'/Annotation_cluster_umap_HPCA.png'),add_lable_plot2,width=12,height=10)
    out_plot_single = SingleR.PlotTsne(Cellcluster_Anno$singler[[1]]$SingleR.single,each_rds@reductions$umap@cell.embeddings, do.labels = F,do.letters = T)
    ggsave(filename=paste0(outdir,'/Single_cell_annot_umap_HPCA.png'),out_plot_single$p,dpi=400,width=10,height=8)
    #ggsave(filename=paste0(outdir,'/Single_cell_annot_umap_HPCA.pdf'),out_plot_single$p,dpi=400,width=10,height=8)
}
#saveRDS(each_rds,file = paste0(outdir,'/singleR_annot1.rds'))
###singler2
clusterAnnotaion=Cellcluster_Anno$singler[[2]]$SingleR.clusters.main$labels #
clusterAnn_df<-data.frame(Annotation=clusterAnnotaion,stringsAsFactors = F)
clusterAnn_df$ID<-rownames(clusterAnn_df) ###
row_length<-nrow(clusterAnnotaion)-1
rownames(clusterAnnotaion)<-paste('cluster',0:row_length,sep='')
clusterAnnotaion<-data.frame(ClusterID=rownames(clusterAnnotaion),Annotation=clusterAnnotaion[,1])
write.table(clusterAnnotaion,file=paste0(outdir,"/Cluster_Cells_Annotation_BLUEENCODE.txt"),quote=F,sep="\t",col.names=T,row.names=F) #save annotation of cluster
singlecellanno<-data.frame(Cell_ID=rownames(Cellcluster_Anno$singler[[2]]$SingleR.single.main$labels),Annotation=Cellcluster_Anno$singler[[2]]$SingleR.single.main$labels[,1], stringsAsFactors = F)##extract cell annotation
singlecellanno$cluster<-mapvalues(singlecellanno$Cell_ID,rownames(each_rds@meta.data),as.character(each_rds@meta.data$seurat_clusters))
write.table(singlecellanno,file=paste0(outdir,"/SinglecellAnnotation_BLUEENCODE.txt"),quote=F,sep="\t",col.names=F,row.names=F) #save cell annotation
###match
each_rds@meta.data$celltype2<-clusterAnn_df$Annotation[match(each_rds@meta.data$seurat_clusters,clusterAnn_df$ID)]
each_rds@meta.data$cluster_celltype2<-paste(each_rds@meta.data$seurat_clusters,each_rds@meta.data$celltype2,sep=':')
each_rds@meta.data$cluster_celltype2<-factor(each_rds@meta.data$cluster_celltype2,levels=unique(each_rds@meta.data[order(each_rds@meta.data$seurat_clusters),]$cluster_celltype2))
###cluster umap/tsne
###single cell umap/tsne
if(cluster=='tsne'){
    add_lable_plot1<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('celltype2'),label = T,raster = F)
    add_lable_plot2<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('cluster_celltype2'),label = T,raster = F)
    ggsave(filename=paste0(outdir,'/Annotation_tsne_BLUEENCODE.png'),add_lable_plot1,width=12,height=10)
    ggsave(filename=paste0(outdir,'/Annotation_cluster_tsne_BLUEENCODE.png'),add_lable_plot2,width=12,height=10)
    out_plot_single = SingleR.PlotTsne(Cellcluster_Anno$singler[[2]]$SingleR.single,each_rds@reductions$tsne@cell.embeddings, do.labels = F,do.letters = T)
    ggsave(filename=paste0(outdir,'/Single_cell_annot_tsne_BLUEENCODE.png'),out_plot_single$p,dpi=400,width=10,height=8)
    #ggsave(filename=paste0(outdir,'/Single_cell_annot_tsne_BLUEENCODE.pdf'),out_plot_single$p,dpi=400,width=10,height=8)
}
if(cluster=='umap'){
    add_lable_plot1<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('celltype2'),label = T,raster = F)
    add_lable_plot2<-DimPlot(each_rds, reduction = cluster,cols = color_panel,group.by = c('cluster_celltype2'),label = T,raster = F)
    ggsave(filename=paste0(outdir,'/Annotation_umap_BLUEENCODE.png'),add_lable_plot1,width=12,height=10)
    ggsave(filename=paste0(outdir,'/Annotation_cluster_umap_BLUEENCODE.png'),add_lable_plot2,width=12,height=10)
    out_plot_single = SingleR.PlotTsne(Cellcluster_Anno$singler[[2]]$SingleR.single,each_rds@reductions$umap@cell.embeddings, do.labels = F,do.letters = T)
    ggsave(filename=paste0(outdir,'/Single_cell_annot_umap_BLUEENCODE.png'),out_plot_single$p,dpi=400,width=10,height=8)
    #ggsave(filename=paste0(outdir,'/Single_cell_annot_umap_BLUEENCODE.pdf'),out_plot_single$p,dpi=400,width=10,height=8)
}
###
#saveRDS(each_rds,file = './singleR_annot2.rds')
