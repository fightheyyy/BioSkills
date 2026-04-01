# .libPaths(c('/share/apps/R/R-4.1.2-4/lib64/R/library'))  # 已注释，使用系统默认R库路径
suppressMessages({
library(dplyr)
library(cowplot)
library(ggplot2)
library(plyr)
library(Seurat)
# library(SeuratWrappers)
library(harmony)
library(stringr)
library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds_files", help="comma-separated RDS file paths")
argv <- add_argument(argv,"--sample_names", help="comma-separated sample names")
argv <- add_argument(argv,"--groups", help="comma-separated group names")
argv <- add_argument(argv,"--outdir", help="output directory")
argv <- add_argument(argv,"--resolution", help="clustering resolution", default=0.5)
argv <- add_argument(argv,"--dims", help="PCA dimensions", default=20)
argv <- parse_args(argv)

rds_files <- strsplit(argv$rds_files, ",")[[1]]
sample_names <- strsplit(argv$sample_names, ",")[[1]]
groups <- strsplit(argv$groups, ",")[[1]]
outdir <- argv$outdir
resolution <- as.numeric(argv$resolution)
dims <- as.numeric(argv$dims)

if(length(rds_files) != length(sample_names) || length(rds_files) != length(groups)){
  stop("rds_files, sample_names and groups must have same length")
}

dir.create(paste0(outdir,'/1.cluster'), recursive=TRUE, showWarnings=FALSE)
dir.create(paste0(outdir,'/merge_rds'), recursive=TRUE, showWarnings=FALSE)

sample_list <- list()
for(i in 1:length(rds_files)){
  sample_list[[i]] <- readRDS(file=rds_files[i])
}

merge_data <- merge(x=sample_list[[1]], y=sample_list[-1], add.cell.ids=sample_names)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merge_data <- merge_data %>% Seurat::NormalizeData(verbose=FALSE) %>% FindVariableFeatures(selection.method='vst', nfeatures=2000) %>% CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes) %>% ScaleData(verbose=FALSE, vars.to.regress=c('S.Score','G2M.Score')) %>% RunPCA(verbose=FALSE)
integrated_data <- merge_data %>% RunHarmony('orig.ident')

sample_group_data <- data_frame(sample=sample_names, group=groups)
integrated_data@meta.data$group <- mapvalues(integrated_data@meta.data$orig.ident, sample_group_data$sample, sample_group_data$group)

integrated_data <- FindNeighbors(integrated_data, reduction='harmony', dims=1:dims)
integrated_data <- FindClusters(integrated_data, resolution=resolution)
integrated_data <- RunTSNE(integrated_data, reduction='harmony', dims=1:dims)
integrated_data <- RunUMAP(integrated_data, reduction='harmony', dims=1:dims)
integrated_data@meta.data$orig.ident <- factor(integrated_data@meta.data$orig.ident, levels=sample_names)
integrated_data@meta.data$group <- factor(integrated_data@meta.data$group, levels=unique(groups))
saveRDS(integrated_data, file=paste0(outdir,'/merge_rds/merge_data.rds'))

color_panel <- c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72','#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B','#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820','#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')

p1 <- DimPlot(integrated_data, label=T, cols=color_panel, raster=F, pt.size=0.4, reduction='umap')
p2 <- DimPlot(integrated_data, label=T, split.by='group', cols=color_panel, raster=F, pt.size=0.4, reduction='umap')
p3 <- DimPlot(integrated_data, label=T, split.by='orig.ident', cols=color_panel, raster=F, pt.size=0.4, ncol=3, reduction='umap')
Idents(integrated_data) <- integrated_data@meta.data$orig.ident
p4 <- DimPlot(integrated_data, pt.size=0.4, raster=F, reduction='umap')
Idents(integrated_data) <- integrated_data@meta.data$group
p44 <- DimPlot(integrated_data, pt.size=0.4, raster=F, reduction='umap')

if("percent.mt" %in% colnames(integrated_data@meta.data)){
  integrated_data[['percent.rp']] <- PercentageFeatureSet(integrated_data, pattern='^RPL|^RPS')
  p5 <- FeaturePlot(integrated_data, features=c('nFeature_RNA','nCount_RNA','percent.mt','percent.rp'), reduction='umap')
}else{
  integrated_data[['percent.mt']] <- PercentageFeatureSet(integrated_data, pattern='^MT-')
  integrated_data[['percent.rp']] <- PercentageFeatureSet(integrated_data, pattern='^RPL|^RPS')
  p5 <- FeaturePlot(integrated_data, features=c('nFeature_RNA','nCount_RNA','percent.mt','percent.rp'), reduction='umap')
}

ggsave(p1, filename=paste0(outdir,'/1.cluster/merge_umap.png'), width=10, height=8)
ggsave(p1, filename=paste0(outdir,'/1.cluster/merge_umap.pdf'), width=10, height=8)
ggsave(p2, filename=paste0(outdir,'/1.cluster/merge_group_umap.png'), width=15, height=12)
ggsave(p2, filename=paste0(outdir,'/1.cluster/merge_group_umap.pdf'), width=15, height=12)
ggsave(p3, filename=paste0(outdir,'/1.cluster/merge_sample_umap.png'), width=15, height=12)
ggsave(p3, filename=paste0(outdir,'/1.cluster/merge_sample_umap.pdf'), width=15, height=12)
ggsave(p4, filename=paste0(outdir,'/1.cluster/sample_umap.png'), width=10, height=8)
ggsave(p4, filename=paste0(outdir,'/1.cluster/sample_umap.pdf'), width=10, height=8)
ggsave(p44, filename=paste0(outdir,'/1.cluster/group_umap.png'), width=10, height=8)
ggsave(p44, filename=paste0(outdir,'/1.cluster/group_umap.pdf'), width=10, height=8)
ggsave(p5, filename=paste0(outdir,'/1.cluster/nFeature_nCount_mt_featureplot.png'), width=10, height=8)
ggsave(p5, filename=paste0(outdir,'/1.cluster/nFeature_nCount_mt_featureplot.pdf'), width=10, height=8)

Idents(integrated_data) <- integrated_data@meta.data$seurat_clusters
p6 <- VlnPlot(integrated_data, features=c('nFeature_RNA','nCount_RNA','percent.mt'), pt.size=0, ncol=1)
ggsave(p6, filename=paste0(outdir,'/1.cluster/nFeature_nCount_mt_violin.png'), width=10, height=12)
ggsave(p6, filename=paste0(outdir,'/1.cluster/nFeature_nCount_mt_violin.pdf'), width=10, height=12)
