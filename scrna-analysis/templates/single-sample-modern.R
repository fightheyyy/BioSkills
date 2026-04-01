# 单细胞自动注释 - 基于现代SingleR API
suppressMessages({
library(dplyr)
library(ggplot2)
library(Seurat)
library(SingleR)
library(celldex)
library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="Seurat RDS文件")
argv <- add_argument(argv,"--species", help="物种: Human/Mouse")
argv <- add_argument(argv,"--reference", help="参考数据集 (可选)", default="auto")
argv <- add_argument(argv,"--outdir", help="输出目录")
argv <- parse_args(argv)

rds_file <- argv$rds
species <- argv$species
ref_name <- argv$reference
outdir <- argv$outdir

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cat("读取Seurat对象...\n")
seurat_obj <- readRDS(rds_file)

# Seurat v5兼容
if(class(seurat_obj@assays$RNA)[1] == "Assay5"){
  seurat_obj <- JoinLayers(seurat_obj)
}

# 提取归一化数据
norm_data <- GetAssayData(seurat_obj, assay="RNA", layer="data")

# 加载参考数据集（使用两个参考）
cat("加载参考数据集...\n")
if(species == "Human"){
  ref1 <- celldex::HumanPrimaryCellAtlasData()
  ref2 <- celldex::BlueprintEncodeData()
}else if(species == "Mouse"){
  ref1 <- celldex::MouseRNAseqData()
  ref2 <- celldex::ImmGenData()
}else{
  stop("物种必须是 Human 或 Mouse")
}

# SingleR注释 - 使用两个参考数据集
cat("执行SingleR注释（参考1）...\n")
# 参考1: Cluster级别
cluster_expr <- AverageExpression(seurat_obj, group.by="seurat_clusters", assays="RNA", layer="data")$RNA
pred_cluster1 <- SingleR(test=cluster_expr, ref=ref1, labels=ref1$label.main)
pred_single1 <- SingleR(test=norm_data, ref=ref1, labels=ref1$label.main)

cat("执行SingleR注释（参考2）...\n")
# 参考2: Cluster级别
pred_cluster2 <- SingleR(test=cluster_expr, ref=ref2, labels=ref2$label.main)
pred_single2 <- SingleR(test=norm_data, ref=ref2, labels=ref2$label.main)

# 处理参考1的结果
cluster_names <- gsub("^g", "", rownames(pred_cluster1))
cluster_map1 <- data.frame(cluster=cluster_names, celltype=pred_cluster1$labels)
celltype_vec1 <- cluster_map1$celltype[match(as.character(seurat_obj$seurat_clusters), cluster_map1$cluster)]
seurat_obj@meta.data$celltype1 <- celltype_vec1
seurat_obj@meta.data$celltype_single1 <- pred_single1$labels

# 处理参考2的结果
cluster_map2 <- data.frame(cluster=cluster_names, celltype=pred_cluster2$labels)
celltype_vec2 <- cluster_map2$celltype[match(as.character(seurat_obj$seurat_clusters), cluster_map2$cluster)]
seurat_obj@meta.data$celltype2 <- celltype_vec2
seurat_obj@meta.data$celltype_single2 <- pred_single2$labels

# 保存注释结果
write.table(cluster_map1, file=paste0(outdir,"/Cluster_Cells_Annotation_HPCA.txt"),
            sep="\t", quote=F, row.names=F, col.names=T)
write.table(cluster_map2, file=paste0(outdir,"/Cluster_Cells_Annotation_BLUEENCODE.txt"),
            sep="\t", quote=F, row.names=F, col.names=T)

single_anno1 <- data.frame(Cell_ID=colnames(seurat_obj), Annotation=pred_single1$labels, cluster=seurat_obj$seurat_clusters)
write.table(single_anno1, file=paste0(outdir,"/SinglecellAnnotation_HPCA.txt"),
            sep="\t", quote=F, row.names=F, col.names=F)

single_anno2 <- data.frame(Cell_ID=colnames(seurat_obj), Annotation=pred_single2$labels, cluster=seurat_obj$seurat_clusters)
write.table(single_anno2, file=paste0(outdir,"/SinglecellAnnotation_BLUEENCODE.txt"),
            sep="\t", quote=F, row.names=F, col.names=F)

# 创建cluster:celltype标签（使用参考1）
seurat_obj$cluster_celltype1 <- paste0(seurat_obj$seurat_clusters, ":", seurat_obj$celltype1)

# 可视化（使用参考1的结果）
color_panel <- c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
                 '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
                 '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
                 '#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')

p1 <- DimPlot(seurat_obj, group.by="celltype1", label=TRUE, cols=color_panel, raster=FALSE)
ggsave(paste0(outdir,"/Annotation_umap_HPCA.png"), p1, width=12, height=10)

p2 <- DimPlot(seurat_obj, group.by="cluster_celltype1", label=TRUE, cols=color_panel, raster=FALSE)
ggsave(paste0(outdir,"/Annotation_cluster_umap_HPCA.png"), p2, width=12, height=10)

p3 <- DimPlot(seurat_obj, group.by="celltype_single1", label=FALSE, cols=color_panel, raster=FALSE)
ggsave(paste0(outdir,"/Single_cell_annot_umap_HPCA.png"), p3, width=10, height=8, dpi=400)

p4 <- DimPlot(seurat_obj, group.by="celltype2", label=TRUE, cols=color_panel, raster=FALSE)
ggsave(paste0(outdir,"/Annotation_umap_BLUEENCODE.png"), p4, width=12, height=10)

p5 <- DimPlot(seurat_obj, group.by="cluster_celltype2", label=TRUE, cols=color_panel, raster=FALSE)
ggsave(paste0(outdir,"/Annotation_cluster_umap_BLUEENCODE.png"), p5, width=12, height=10)

p6 <- DimPlot(seurat_obj, group.by="celltype_single2", label=FALSE, cols=color_panel, raster=FALSE)
ggsave(paste0(outdir,"/Single_cell_annot_umap_BLUEENCODE.png"), p6, width=10, height=8, dpi=400)

cat("注释完成！\n")
