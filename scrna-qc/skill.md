---
name: scrna-qc
description: "单细胞RNA-seq质控过滤。支持多样本批量处理，可调整质控参数（min/max features, 线粒体比例）。"
category: 科研
invocable: user
autoInvocable: true
argument-hint: "<input-dirs> <sample-names>"
---

# 单细胞质控过滤

对10X Genomics原始数据进行质控过滤。

## 功能

- 读取10X filtered_feature_bc_matrix数据
- 过滤低质量细胞（基于基因数、UMI数、线粒体比例）
- 支持Human/Mouse/Rat物种
- 输出过滤后的RDS文件

## 参数

```bash
Rscript filter-qc.R \
  --input_dirs "dir1,dir2,dir3" \
  --sample_names "s1,s2,s3" \
  --species Human \
  --outdir /path/out \
  --min_features 400 \
  --max_features 6000 \
  --max_mt 25
```

**必需参数：**
- `--input_dirs`: 10X数据目录（逗号分隔）
- `--sample_names`: 样本名称（逗号分隔）
- `--species`: 物种（Human/Mouse/Rat）
- `--outdir`: 输出目录

**可选参数：**
- `--min_features`: 最小基因数（默认400）
- `--max_features`: 最大基因数（默认6000）
- `--max_mt`: 最大线粒体比例（默认25）

## 参数调整建议

**min_features/max_features:**
- 低质量数据：min=200, max=4000
- 高质量数据：min=500, max=8000
- 查看小提琴图后调整

**max_mt:**
- 组织样本：15-25%
- 细胞系：5-10%
- 应激/凋亡样本：可放宽至30%

## 输出

- `sample_rds/`: 过滤后的RDS文件
- `qc_plots/`: 质控图（小提琴图、散点图）

## 依赖

```r
install.packages(c("dplyr", "ggplot2", "cowplot", "patchwork", "argparser", "stringr"))
BiocManager::install("Seurat")
```
