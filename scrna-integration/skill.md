---
name: scrna-integration
description: "单细胞样本整合与聚类分析。包含Harmony批次校正、Marker基因鉴定、细胞数量统计。"
category: 科研
invocable: user
autoInvocable: true
argument-hint: "<rds-files> <sample-names>"
---

# 单细胞整合与聚类分析

多样本整合、Harmony批次校正、聚类、Marker基因鉴定、细胞统计。

## 包含的分析

1. **merge-harmony.R** - 样本合并与Harmony批次校正
2. **marker-gene.R** - Marker基因鉴定
3. **cell-number.R** - 细胞数量统计

## 使用流程

### 1. 样本合并与批次校正

```bash
Rscript merge-harmony.R \
  --rds_files "s1.rds,s2.rds" \
  --sample_names "s1,s2" \
  --groups "ctrl,treat" \
  --outdir /path/out \
  --resolution 0.5 \
  --dims 20
```

**参数：**
- `--rds_files`: RDS文件路径（逗号分隔）
- `--sample_names`: 样本名称
- `--groups`: 分组信息
- `--outdir`: 输出目录
- `--resolution`: 聚类分辨率（默认0.5，范围0.1-2.0）
- `--dims`: PCA维度（默认20）

**输出：**
- `merge_rds/merge_data.rds`: 合并后的数据
- `1.cluster/`: UMAP/tSNE图

### 2. Marker基因鉴定

**增强版脚本**: `marker-gene-enhanced.R`

```bash
Rscript marker-gene-enhanced.R \
  --rds merge_data.rds \
  --outdir /path/out \
  --reduction umap \
  --species Mouse \
  --top_n 10 \
  --custom_genes "Cd3d,Cd8a,Cd4,Foxp3" \
  --colors "#FF0000,#00FF00,#0000FF"
```

**新增参数：**
- `--top_n`: 每个cluster的top N marker基因（默认10）
- `--custom_genes`: 自定义基因列表（逗号分隔）
- `--colors`: 自定义颜色（逗号分隔）

**输出：**
- `2.marker/all_markgene.xls`: 所有Marker基因
- `2.marker/sig_markgene.xls`: 显著Marker基因
- `2.marker/topN_markers.xls`: Top N marker基因表
- `2.marker/heatmap/`: 热图
- `2.marker/violin/`: 小提琴图（无散点）
- `2.marker/featureplot/`: 特征图（4列排列）
- `2.marker/custom_genes/`: 自定义基因可视化

### 3. 细胞数量统计

```bash
Rscript cell-number.R \
  --rds1 merge_data.rds \
  --outdir /path/out \
  --group group \
  --celltype seurat_clusters
```

**输出：**
- `4.cell_number/`: 细胞数量统计表和柱状图

## 参数调整建议

**resolution（聚类分辨率）：**
- 0.1-0.5: 粗聚类，适合大类群鉴定
- 0.5-1.0: 标准聚类
- 1.0-2.0: 细聚类，适合亚群分析

**dims（PCA维度）：**
- 查看ElbowPlot确定
- 通常15-30维

## 常见需求

**需求1**: 筛选top N marker基因
```bash
--top_n 20  # 每个cluster的top 20基因
```

**需求2**: Marker基因小提琴图和FeaturePlot
- 自动生成，小提琴图无散点，FeaturePlot按4列排列

**需求3**: 自定义分群颜色
```bash
--colors "#FF0000,#00FF00,#0000FF,#FFFF00"
```

**需求11**: 自定义基因列表可视化
```bash
--custom_genes "Cd3d,Cd8a,Cd4,Foxp3,Il2ra,Ctla4"
```
输出到 `custom_genes/` 目录

## 依赖

```r
install.packages(c("dplyr", "ggplot2", "cowplot", "patchwork", "argparser", "stringr", "plyr"))
BiocManager::install(c("Seurat", "harmony"))
```
