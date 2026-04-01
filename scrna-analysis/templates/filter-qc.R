# .libPaths(c('/share/apps/R/R-4.1.2-4/lib64/R/library'))  # 已注释，使用系统默认R库路径
suppressMessages({
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(Seurat)
library(stringr)
library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--input_dirs", help="comma-separated 10X data directories")
argv <- add_argument(argv,"--sample_names", help="comma-separated sample names")
argv <- add_argument(argv,"--species", help="species: Human/Mouse/Rat")
argv <- add_argument(argv,"--outdir", help="output directory")
argv <- add_argument(argv,"--min_features", help="min nFeature_RNA", default=400)
argv <- add_argument(argv,"--max_features", help="max nFeature_RNA", default=6000)
argv <- add_argument(argv,"--max_mt", help="max percent.mt", default=25)
argv <- parse_args(argv)

input_dirs <- strsplit(argv$input_dirs, ",")[[1]]
sample_names <- strsplit(argv$sample_names, ",")[[1]]
species <- argv$species
outdir <- argv$outdir
min_features <- as.numeric(argv$min_features)
max_features <- as.numeric(argv$max_features)
max_mt <- as.numeric(argv$max_mt)

if(length(input_dirs) != length(sample_names)){
  stop("input_dirs and sample_names must have same length")
}

# 根据物种设置线粒体基因模式
if(species == 'Human'){
  mt_pattern <- '^MT-'
}else if(species == 'Mouse'){
  mt_pattern <- '^mt-'
}else if(species == 'Rat'){
  mt_pattern <- '^Mt-'
}else{
  stop("species must be Human, Mouse, or Rat")
}

dir.create(paste0(outdir,'/0.QC'), recursive=TRUE, showWarnings=FALSE)
dir.create(paste0(outdir,'/sample_rds'), recursive=TRUE, showWarnings=FALSE)

sample_list <- list()
nFeature_plots <- list()
nCount_plots <- list()
mt_plots <- list()

for(i in 1:length(input_dirs)){
  sample_data <- Read10X(input_dirs[i])
  sample_obj <- CreateSeuratObject(counts=sample_data, project=sample_names[i], min.cells=3, min.features=200)
  sample_obj[['percent.mt']] <- PercentageFeatureSet(sample_obj, pattern=mt_pattern)
  sample_obj[['percent.rp']] <- PercentageFeatureSet(sample_obj, pattern='^RPL|^RPS')
  nFeature_plots[[i]] <- VlnPlot(sample_obj, features='nFeature_RNA', pt.size=0)
  nCount_plots[[i]] <- VlnPlot(sample_obj, features='nCount_RNA', pt.size=0)
  mt_plots[[i]] <- VlnPlot(sample_obj, features='percent.mt', pt.size=0)
  sample_list[[i]] <- sample_obj
}

plot_nFeature <- wrap_plots(nFeature_plots, ncol=3)
plot_nCount <- wrap_plots(nCount_plots, ncol=3)
plot_mt <- wrap_plots(mt_plots, ncol=3)
ggsave(plot_nFeature, filename=paste0(outdir,'/0.QC/nFeature_raw.png'), width=20, height=15)
ggsave(plot_nCount, filename=paste0(outdir,'/0.QC/nCount_raw.png'), width=20, height=15)
ggsave(plot_mt, filename=paste0(outdir,'/0.QC/mt_raw.png'), width=20, height=15)

nFeature_plots <- list()
nCount_plots <- list()
mt_plots <- list()

for(i in 1:length(sample_list)){
  print(paste0(sample_names[i], " raw cell number:"))
  print(ncol(sample_list[[i]]))
  sample_list[[i]] <- subset(sample_list[[i]], subset=nFeature_RNA>min_features & nFeature_RNA<max_features & percent.mt<max_mt)
  print(paste0(sample_names[i], " filter cell number:"))
  print(ncol(sample_list[[i]]))
  saveRDS(sample_list[[i]], file=paste0(outdir,'/sample_rds/',sample_names[i],'.rds'))
  nFeature_plots[[i]] <- VlnPlot(sample_list[[i]], features='nFeature_RNA', pt.size=0)
  nCount_plots[[i]] <- VlnPlot(sample_list[[i]], features='nCount_RNA', pt.size=0)
  mt_plots[[i]] <- VlnPlot(sample_list[[i]], features='percent.mt', pt.size=0)
}

plot_nFeature <- wrap_plots(nFeature_plots, ncol=3)
plot_nCount <- wrap_plots(nCount_plots, ncol=3)
plot_mt <- wrap_plots(mt_plots, ncol=3)
ggsave(plot_nFeature, filename=paste0(outdir,'/0.QC/nFeature_filter.png'), width=20, height=15)
ggsave(plot_nCount, filename=paste0(outdir,'/0.QC/nCount_filter.png'), width=20, height=15)
ggsave(plot_mt, filename=paste0(outdir,'/0.QC/mt_filter.png'), width=20, height=15)
