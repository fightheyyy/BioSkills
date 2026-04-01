# 细胞比例可视化脚本
suppressMessages({
library(dplyr)
library(ggplot2)
library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="Seurat RDS文件")
argv <- add_argument(argv,"--outdir", help="输出目录")
argv <- add_argument(argv,"--celltype", help="细胞类型列名", default="seurat_clusters")
argv <- add_argument(argv,"--group", help="分组列名(可选)", default=NA)
argv <- parse_args(argv)

rds_file <- argv$rds
outdir <- argv$outdir
celltype_col <- argv$celltype
group_col <- argv$group

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cat("读取Seurat对象...\n")
seurat_obj <- readRDS(rds_file)

# 计算细胞类型比例
celltype_counts <- table(seurat_obj@meta.data[[celltype_col]])
celltype_prop <- data.frame(
  CellType = names(celltype_counts),
  Count = as.numeric(celltype_counts),
  Proportion = as.numeric(celltype_counts) / sum(celltype_counts)
)
celltype_prop <- celltype_prop[order(celltype_prop$Count, decreasing=TRUE),]
celltype_prop$CellType <- factor(celltype_prop$CellType, levels=celltype_prop$CellType)

write.table(celltype_prop, paste0(outdir,'/celltype_proportion.txt'), sep="\t", quote=F, row.names=F)

# 需求9: 细胞类型占比柱状图（从大到小排序）
cat("生成细胞类型占比柱状图...\n")
p1 <- ggplot(celltype_prop, aes(x=CellType, y=Proportion, fill=CellType)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
  labs(x="Cell Type", y="Proportion", title="Cell Type Proportion") +
  scale_y_continuous(labels=scales::percent)
ggsave(p1, filename=paste0(outdir,'/celltype_proportion_barplot.png'), width=10, height=6)

# 需求10: 分组堆叠图
if(!is.na(group_col) && group_col %in% colnames(seurat_obj@meta.data)){
  cat("生成分组堆叠图...\n")
  group_celltype <- table(seurat_obj@meta.data[[group_col]], seurat_obj@meta.data[[celltype_col]])
  group_celltype_df <- as.data.frame(group_celltype)
  colnames(group_celltype_df) <- c("Group", "CellType", "Count")

  group_celltype_df <- group_celltype_df %>%
    group_by(Group) %>%
    mutate(Proportion = Count / sum(Count))

  write.table(group_celltype_df, paste0(outdir,'/group_celltype_proportion.txt'), sep="\t", quote=F, row.names=F)

  p2 <- ggplot(group_celltype_df, aes(x=Group, y=Proportion, fill=CellType)) +
    geom_bar(stat="identity", position="stack") +
    theme_classic() +
    labs(x="Group", y="Proportion", title="Cell Type Proportion by Group") +
    scale_y_continuous(labels=scales::percent)
  ggsave(p2, filename=paste0(outdir,'/group_celltype_stacked.png'), width=10, height=6)
}

cat("细胞比例分析完成！\n")

