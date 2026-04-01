# 参数化模板使用说明

## 全部 6 个参数化脚本

### 1. seurat_merge_cluster_harmony.R
多样本合并与 Harmony 批次校正

**参数：** 需要修改脚本中的样本路径、样本名、分组信息

### 2. markgene.R
Marker 基因鉴定

**用法：**
```bash
Rscript markgene.R \
  --rds "merged_data.rds" \
  --outdir "/path/to/output" \
  --reduction "umap" \
  --species "Human"
```

### 3. diff_enrich_volcano_4.1.R
差异分析与富集分析

**用法：**
```bash
Rscript diff_enrich_volcano_4.1.R \
  --rds "data.rds" \
  --group "AvsB" \
  --group_flag "group" \
  --outdir "/path/to/output" \
  --species "Human" \
  --volcano_flag 1 \
  --p_val_adj 0.05 \
  --avg_log2FC 0.25
```

### 4. seurat_cell_numer.R
细胞数量统计

**用法：**
```bash
Rscript seurat_cell_numer.R \
  --rds1 "merged_data.rds" \
  --group "group" \
  --celltype "seurat_clusters" \
  --outdir "/path/to/output"
```

### 5. seurat_filter_human.R
质控过滤

**参数：** 需要修改脚本中的数据路径和过滤阈值

### 6. single.R
单样本完整分析

**用法：**
```bash
Rscript single.R \
  --rds "sample.rds" \
  --cluster "umap" \
  --speices "Human" \
  --outdir "/path/to/output"
```
