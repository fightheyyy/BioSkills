---
name: scrna-diff
description: "单细胞差异表达分析与GO富集。支持组间比较、火山图、GO富集分析（KEGG已禁用）。"
invocable: user
autoInvocable: true
argument-hint: "<rds-file> <group1vsgroup2>"
---

# 单细胞差异分析与富集

组间差异表达分析与GO功能富集。

## 功能

- 组间差异基因鉴定
- 火山图可视化
- GO富集分析（BP/CC/MF）
- 支持上调/下调基因分别富集

## 使用方法

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

## 依赖

```r
install.packages(c("dplyr", "ggplot2", "cowplot", "argparser", "stringr", "ggrepel"))
BiocManager::install(c("Seurat", "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db"))
```
