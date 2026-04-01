---
name: scrna-integration
description: "单细胞样本整合与聚类分析。包含Harmony批次校正、Marker基因鉴定、细胞数量统计。"
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

```bash
Rscript marker-gene.R \
  --rds merge_data.rds \
  --outdir /path/out \
  --reduction umap \
  --species Mouse
```

**输出：**
- `2.marker/all_markgene.xls`: 所有Marker基因
- `2.marker/sig_markgene.xls`: 显著Marker基因
- `2.marker/heatmap/`: 热图
- `2.marker/violin/`: 小提琴图
- `2.marker/featureplot/`: 特征图

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

## 依赖

```r
install.packages(c("dplyr", "ggplot2", "cowplot", "patchwork", "argparser", "stringr", "plyr"))
BiocManager::install(c("Seurat", "harmony"))
```
