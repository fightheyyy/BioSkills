---
name: scrna-diff
description: "单细胞差异表达分析与GO富集。支持组间比较、火山图、GO富集分析、细胞比例可视化。"
category: 科研
invocable: user
autoInvocable: true
argument-hint: "<rds-file> <group1vsgroup2>"
---

# 单细胞差异分析与富集

组间差异表达分析、GO功能富集、细胞比例统计与可视化。

## 包含的分析

1. **diff-enrich.R** - 差异分析与GO富集
2. **cell-proportion.R** - 细胞比例可视化（新增）

## 使用方法

### 1. 差异分析与GO富集

```bash
Rscript diff-enrich.R \
  --rds merge_data.rds \
  --group "group1vsgroup2" \
  --group_flag "group" \
  --species Mouse \
  --outdir /path/out \
  --volcano_flag 0 \
  --p_val_adj 0.05 \
  --avg_log2FC 0.25
```

### 2. 细胞比例可视化（新增）

```bash
Rscript cell-proportion.R \
  --rds merge_data.rds \
  --outdir /path/out \
  --celltype seurat_clusters \
  --group group
```

**参数：**
- `--celltype`: 细胞类型列名（默认seurat_clusters）
- `--group`: 分组列名（可选，用于分组堆叠图）

**输出：**
- `celltype_proportion.txt`: 细胞类型比例表
- `celltype_proportion_barplot.png`: 细胞类型占比柱状图（从大到小排序）
- `group_celltype_stacked.png`: 分组堆叠图（如果提供group参数）

## 参数

**必需：**
- `--rds`: Seurat RDS文件
- `--group`: 比较组（格式：AvsB）
- `--group_flag`: 分组列名（如"group"）
- `--species`: 物种（Human/Mouse/Rat）
- `--outdir`: 输出目录

**可选：**
- `--volcano_flag`: 火山图标注（0=部分基因，1=全部基因）
- `--p_val_adj`: p值阈值（默认0.05）
- `--avg_log2FC`: logFC阈值（默认0.25）

## 参数调整建议

**p_val_adj:**
- 严格筛选：0.01
- 标准：0.05
- 宽松：0.1

**avg_log2FC:**
- 严格：0.5-1.0
- 标准：0.25
- 宽松：0.1

## 输出

- `diff/`: 差异基因表格
- `diff/volcano/`: 火山图
- `enrich/ALL/`: 全部差异基因GO富集
- `enrich/UP/`: 上调基因GO富集
- `enrich/DOWN/`: 下调基因GO富集

## 注意

- KEGG富集已禁用（KEGG.db在Bioconductor 3.22+不可用）
- 仅提供GO富集分析

## 常见需求

**需求4**: 特定cluster的富集分析
- 先用diff-enrich.R找到差异基因，再对特定cluster进行富集

**需求9**: 细胞类型占比柱状图
```bash
Rscript cell-proportion.R --rds data.rds --outdir /path/out --celltype seurat_clusters
```
自动按从大到小排序

**需求10**: 分组堆叠图
```bash
Rscript cell-proportion.R --rds data.rds --outdir /path/out --celltype seurat_clusters --group group
```

## 依赖

```r
install.packages(c("dplyr", "ggplot2", "cowplot", "argparser", "stringr", "ggrepel"))
BiocManager::install(c("Seurat", "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db"))
```
