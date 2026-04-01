# 增强版Marker基因分析脚本
suppressMessages({
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
library(argparser)})
plan("multicore", workers = 8)
options(future.globals.maxSize = 1000000 * 1024^5)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="Seurat RDS文件")
argv <- add_argument(argv,"--outdir", help="输出目录")
argv <- add_argument(argv,"--reduction", help="降维方法: umap/tsne", default="umap")
argv <- add_argument(argv,"--species", help="物种: Human/Mouse/Rat")
argv <- add_argument(argv,"--top_n", help="每个cluster的top N marker基因", default=10)
argv <- add_argument(argv,"--custom_genes", help="自定义基因列表(逗号分隔)", default=NA)
argv <- add_argument(argv,"--colors", help="自定义颜色(逗号分隔)", default=NA)
argv <- parse_args(argv)

rds <- argv$rds
outdir <- argv$outdir
reduction <- argv$reduction
species <- argv$species
top_n <- as.numeric(argv$top_n)
custom_genes <- argv$custom_genes
custom_colors <- argv$colors

cat("读取Seurat对象...\n")
integrated_data <- readRDS(rds)
DefaultAssay(integrated_data) <- 'RNA'
integrated_data$seurat_clusters <- factor(integrated_data$seurat_clusters, levels=unique(sort(integrated_data$seurat_clusters)))
Idents(integrated_data) <- integrated_data@meta.data$seurat_clusters
integrated_data <- NormalizeData(integrated_data, verbose=FALSE)
integrated_data <- JoinLayers(integrated_data)

cat("寻找Marker基因...\n")
integrated_data.allmarkers <- FindAllMarkers(integrated_data, only.pos=TRUE, min.pct=0.2, logfc.threshold=0.25)
integrated_data.allmarkers <- integrated_data.allmarkers[!grepl(integrated_data.allmarkers$gene, pattern='^MT-|^mt-|^Mt-'),]

all_sig_marker_gene <- integrated_data.allmarkers[integrated_data.allmarkers$p_val_adj<0.05,]
write.table(integrated_data.allmarkers, paste0(outdir,'/all_markgene.xls'), quote=F, row.names=F, sep='\t')
write.table(all_sig_marker_gene, paste0(outdir,'/sig_markgene.xls'), quote=F, row.names=F, sep='\t')

all.genes <- rownames(integrated_data)
integrated_data <- ScaleData(integrated_data, features=all.genes)

# 提取top N marker基因
topN <- all_sig_marker_gene %>% group_by(cluster) %>% top_n(n=top_n, wt=avg_log2FC)
topN_plus2 <- all_sig_marker_gene %>% group_by(cluster) %>% top_n(n=top_n+2, wt=avg_log2FC)
write.table(topN, paste0(outdir,'/top',top_n,'_markers.xls'), quote=F, row.names=F, sep='\t')

# 创建输出目录
dir.create(paste0(outdir,'/heatmap'), showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(outdir,'/featureplot'), showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(outdir,'/violin'), showWarnings=FALSE, recursive=TRUE)

# 设置颜色
if(!is.na(custom_colors)){
  color_panel <- strsplit(custom_colors, ",")[[1]]
}else{
  color_panel <- c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
                   '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
                   '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
                   '#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')
}

cat("生成热图...\n")
p_heatmap <- DoHeatmap(integrated_data, features=as.character(topN$gene)) +
  scale_fill_gradient2(low="#0000FF", mid="#FFF8DC", high="#FF0000")
ggsave(p_heatmap, filename=paste0(outdir,'/heatmap/top',top_n,'_heatmap.png'), width=25, height=20)

cat("生成各cluster的小提琴图和FeaturePlot...\n")
for(i in levels(integrated_data@meta.data$seurat_clusters)){
  cluster_genes <- as.character(topN_plus2[topN_plus2$cluster==i,]$gene)
  if(length(cluster_genes) > 0){
    p_vln <- VlnPlot(integrated_data, features=cluster_genes, pt.size=0, cols=color_panel)
    ggsave(p_vln, filename=paste0(outdir,'/violin/violin_cluster',i,'.png'), width=20, height=15)

    p_feat <- FeaturePlot(integrated_data, features=cluster_genes, label=TRUE, reduction=reduction, raster=FALSE, ncol=4) &
      scale_colour_gradientn(colours=c("lightblue","orange","red"))
    ggsave(p_feat, filename=paste0(outdir,'/featureplot/featureplot_cluster',i,'.png'), width=20, height=15)
  }
}

# 自定义基因列表可视化
if(!is.na(custom_genes)){
  cat("生成自定义基因列表可视化...\n")
  gene_list <- strsplit(custom_genes, ",")[[1]]
  gene_list <- gene_list[gene_list %in% rownames(integrated_data)]

  if(length(gene_list) > 0){
    dir.create(paste0(outdir,'/custom_genes'), showWarnings=FALSE, recursive=TRUE)
    p_custom_vln <- VlnPlot(integrated_data, features=gene_list, pt.size=0, cols=color_panel, ncol=4)
    ggsave(p_custom_vln, filename=paste0(outdir,'/custom_genes/custom_violin.png'), width=20, height=15)

    p_custom_feat <- FeaturePlot(integrated_data, features=gene_list, label=TRUE, reduction=reduction, raster=FALSE, ncol=4) &
      scale_colour_gradientn(colours=c("lightblue","orange","red"))
    ggsave(p_custom_feat, filename=paste0(outdir,'/custom_genes/custom_featureplot.png'), width=20, height=15)
  }
}

cat("Marker基因分析完成！\n")

