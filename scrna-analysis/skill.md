<!-- ---
name: scrna-analysis
description: "单细胞RNA-seq分析工具。基于Seurat的标准分析流程，支持：多样本合并与Harmony批次校正、Marker基因鉴定、差异分析与富集、细胞数量统计、质控过滤、单样本分析。"
invocable: user
autoInvocable: true
argument-hint: "<analysis-type>"
--- -->

# 单细胞 RNA-seq 分析

基于 Seurat 的单细胞分析标准流程。

## 支持的分析类型

1. **merge-harmony** - 多样本合并与 Harmony 批次校正
2. **marker-gene** - Marker 基因鉴定
3. **diff-enrich** - 差异分析与富集分析
4. **cell-number** - 细胞数量统计
5. **filter-qc** - 质控过滤
6. **single-sample** - 单样本完整分析

## 工作流程

当用户请求单细胞分析时：

1. 询问分析类型和参数
2. 从 templates/ 读取对应的 R 脚本模板
3. 根据用户输入替换参数（使用 sed 或直接生成）
4. 生成可执行的 R 脚本到用户指定目录
5. 询问是否立即执行

## 参数示例

### 合并分析 (merge-harmony)
```
样本路径: /path/to/sample1.rds, /path/to/sample2.rds
样本名称: sample1, sample2
分组信息: control, treatment
输出目录: /path/to/output
分辨率: 0.5 (可选)
PCA维度: 1:20 (可选)
```

### Marker 基因 (marker-gene)
```
输入RDS: /path/to/merged.rds
物种: Human/Mouse/Rat
输出目录: /path/to/output
降维方法: umap/tsne (可选)
```

### 差异分析 (diff-enrich)
```
输入RDS: /path/to/data.rds
比较组: group1 vs group2
logFC阈值: 0.25 (可选)
p值阈值: 0.05 (可选)
输出目录: /path/to/output
```

## 可用脚本

templates/ 目录包含 5 个参数化 R 脚本，均支持命令行参数：

### 1. filter-qc.R - 质控过滤
```bash
Rscript filter-qc.R --input_dirs "dir1,dir2,dir3" --sample_names "s1,s2,s3" --species Human --outdir /path/out --min_features 400 --max_features 6000 --max_mt 25
```
**重要**: 必须指定 `--species` 参数（Human/Mouse/Rat），用于正确识别线粒体基因

### 2. merge-harmony.R - 多样本合并与Harmony批次校正
```bash
Rscript merge-harmony.R --rds_files "s1.rds,s2.rds" --sample_names "s1,s2" --groups "ctrl,treat" --outdir /path/out --resolution 0.5 --dims 20
```

### 3. marker-gene.R - Marker基因鉴定
```bash
Rscript marker-gene.R --rds data.rds --outdir /path/out --reduction umap --species Human --gene_annot_dir /path/to/genome
```

### 4. cell-number.R - 细胞数量统计
```bash
Rscript cell-number.R --rds1 merge.rds --outdir /path/out --group group --celltype seurat_clusters
```

### 5. single-sample.R - 单样本注释分析
```bash
Rscript single-sample.R --rds data.rds --cluster umap --speices Human --outdir /path/out
```

## 依赖环境

### R包依赖
分析前需确保已安装以下R包：

**核心包：**
- **Seurat** - 单细胞分析核心包
- **harmony** - 批次校正
- **argparser** - 命令行参数解析

**可视化包：**
- **dplyr, ggplot2, cowplot, patchwork** - 数据处理与可视化

**可选包（特定分析需要）：**
- **SingleR** - 单样本细胞类型注释（single-sample.R需要）
- **clusterProfiler, org.Mm.eg.db, org.Hs.eg.db** - 富集分析（diff-enrich.R需要）

安装命令：
```r
# 基础包
install.packages(c("dplyr", "ggplot2", "cowplot", "patchwork", "argparser", "stringr", "plyr"))

# Bioconductor包
install.packages("BiocManager")
BiocManager::install(c("Seurat", "harmony"))

# 可选：单样本注释
BiocManager::install("SingleR")

# 可选：富集分析
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db"))
```

## 推荐分析流程

标准单细胞分析的执行顺序：

1. **filter-qc.R** - 质控过滤（第一步）
   - 读取10X原始数据
   - 过滤低质量细胞
   - 输出过滤后的RDS文件

2. **merge-harmony.R** - 多样本合并（第二步）
   - 读取过滤后的RDS文件
   - Harmony批次校正
   - 聚类分析
   - 输出合并后的RDS文件

3. **后续分析**（可并行执行）：
   - **marker-gene.R** - 鉴定各cluster的标记基因
   - **cell-number.R** - 统计细胞数量和比例
   - **diff-enrich.R** - 差异分析与功能富集
   - **single-sample.R** - 单样本细胞类型注释

## 注意事项

- 所有脚本已移除硬编码的R库路径，使用系统默认路径
- marker-gene.R 的基因注释文件路径已参数化（可选参数 --gene_annot_dir）
- filter-qc.R 和 merge-harmony.R 已完全参数化，支持任意数量样本
- 建议按照推荐流程顺序执行，确保输入输出文件路径正确衔接
