# 单细胞 RNA-seq 分析 Skill

## 模板说明

本 skill 包含 6 个 R 脚本模板，覆盖单细胞分析的常见场景：

### 1. seurat_merge_cluster_harmony.R
**功能：** 多样本合并与 Harmony 批次校正
- 合并多个 Seurat 对象
- 使用 Harmony 消除批次效应
- 聚类和降维（PCA、UMAP、tSNE）
- 生成完整可视化结果（cluster图、分组图、样本图、QC图）

### 2. markgene.R
**功能：** Marker 基因鉴定
- FindAllMarkers 找各 cluster 的差异基因
- 过滤线粒体基因
- 添加基因注释（支持 Human/Mouse/Rat）
- 生成热图、小提琴图、FeaturePlot

### 3. diff_enrich_volcano_4.1.R
**功能：** 差异分析与富集分析
- 两组间差异表达分析
- GO/KEGG 富集分析
- 火山图可视化
- 富集结果可视化

### 4. seurat_cell_numer.R
**功能：** 细胞数量统计
- 统计各 cluster 的细胞数
- 按样本分组统计
- 生成统计表格和图表

### 5. seurat_filter_human.R
**功能：** 数据质控与过滤
- 基因数、UMI 数过滤
- 线粒体比例过滤
- 去除双细胞
- 标准化和归一化

### 6. single.R
**功能：** 单样本完整分析流程
- 从原始数据到聚类的完整流程
- 质控、归一化、降维、聚类
- 生成标准分析报告

## 使用方式

XiaoBa 会根据用户需求调用对应脚本，传递参数执行分析。所有脚本输出的图片与原始版本完全一致。
