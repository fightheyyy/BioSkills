# .libPaths('/share/apps/R/R-4.1.2-4/lib64/R/library')  # 已注释，使用系统默认R库路径
suppressMessages({
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
library(argparser)})
plan("multicore", workers = 8)
options(future.globals.maxSize = 1000000 * 1024^5)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="the seurat rds file")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--reduction", help="the reduction of cluster")
argv <- add_argument(argv,"--species", help="the species of project")
argv <- add_argument(argv,"--gene_annot_dir", help="the gene annotation directory (optional)", default=NA)
argv <- parse_args(argv)

rds <- argv$rds
outdir <- argv$outdir
reduction <- argv$reduction
gene_annot_dir <- argv$gene_annot_dir
species <- argv$species


integrated_data<-readRDS(rds)

DefaultAssay(integrated_data)<-'RNA'
#integrated_data$seurat_clusters<-factor(integrated_data$seurat_clusters,levels=unique(integrated_data$seurat_clusters))
#integrated_data$seurat_clusters<-as.character(integrated_data$seurat_clusters)
integrated_data$seurat_clusters<-factor(integrated_data$seurat_clusters,levels=unique(sort(integrated_data$seurat_clusters)))
Idents(integrated_data)<-integrated_data@meta.data$seurat_clusters
integrated_data<-NormalizeData(integrated_data,verbose = FALSE)
integrated_data <- JoinLayers(integrated_data)
integrated_data.allmarkers <- FindAllMarkers(integrated_data, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
integrated_data.allmarkers<-integrated_data.allmarkers[!grepl(integrated_data.allmarkers$gene,pattern = '^MT-|^mt-|^Mt-'),]
if (!is.na(species) && !is.na(gene_annot_dir)){
  if(species == 'Human'){
    gene_annot_file<-read.delim(file='/share/home/liny/Script/bulkseq/genome/Animal/Homo_sapiens/GRCh38/GRCh38_gene.txt',header=T,sep='\t')
  }
  if(species == 'Mouse'){
    gene_annot_file<-read.delim(file='/share/home/liny/Script/bulkseq/genome/Animal/Mus_musculus/Mm10/Mm10_gene.txt',header=T,sep='\t')
  }
  if(species == 'Rat'){
    gene_annot_file<-read.delim(file='/share/home/liny/Script/bulkseq/genome/Animal/Rattus_norvegicus/Rn6/Rn6_gene.txt',header=T,sep='\t')
  }
  gene_annot_file<-gene_annot_file[,c('gene_name','gene_description')]
  gene_annot_file<-unique(gene_annot_file)
  uniq_gene_list<-c()
  uniq_gene_annot_file<-data.frame()
  i=1
  for(each in as.character(gene_annot_file$gene_name)){
    if(!(each %in% uniq_gene_list)){
      uniq_gene_annot_file<-rbind(uniq_gene_annot_file,gene_annot_file[i,])
      uniq_gene_list<-c(uniq_gene_list,each)}
    i<-i+1
  }
  rownames(uniq_gene_annot_file)<-uniq_gene_annot_file$gene_name
  integrated_data.allmarkers$gene_description<-uniq_gene_annot_file[as.character(integrated_data.allmarkers$gene),]$gene_description
}
all_sig_marker_gene<-integrated_data.allmarkers[integrated_data.allmarkers$p_val_adj<0.05,]
write.table(integrated_data.allmarkers,paste(outdir,'/all_markgene.xls',sep = ''),quote = F,row.names = F,sep = '\t')
write.table(all_sig_marker_gene,paste(outdir,'/sig_markgene.xls',sep = ''),quote = F,row.names = F,sep = '\t')
all.genes <- rownames(integrated_data)
integrated_data <- ScaleData(integrated_data, features = all.genes)
top10 <- all_sig_marker_gene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top12 <- all_sig_marker_gene %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC)
#p1<-DoHeatmap(integrated_data, features = as.character(top10$gene)) + NoLegend()+scale_fill_gradient2(low = "#0000FF",mid = "#FFF8DC",high = "#FF0000")
p2<-DoHeatmap(integrated_data, features = as.character(top10$gene)) +scale_fill_gradient2(low = "#0000FF",mid = "#FFF8DC",high = "#FF0000")

if(file.exists(paste0(outdir,'/heatmap'))){
  unlink(paste0(outdir,'/heatmap'), recursive=TRUE)
}else{dir.create(paste0(outdir,'/heatmap'))}
if(file.exists(paste0(outdir,'/featureplot'))){
  unlink(paste0(outdir,'/featureplot'), recursive=TRUE)
}else{dir.create(paste0(outdir,'/featureplot'))}
if(file.exists(paste0(outdir,'/violin'))){
  unlink(paste0(outdir,'/violin'), recursive=TRUE)
}else{dir.create(paste0(outdir,'/violin'))}

dir.create(paste0(outdir,'/heatmap'))
dir.create(paste0(outdir,'/featureplot'))
dir.create(paste0(outdir,'/violin'))

ggsave(p2,filename = paste0(outdir,'/heatmap/top10_heatmap.pdf'),width = 25,height = 20)

for (i in levels(integrated_data@meta.data$seurat_clusters)) {
    p1<-VlnPlot(integrated_data,features = as.character(top12[top12$cluster==i,]$gene),pt.size = 0)
    ggsave(p1,filename = paste0(outdir,'/violin/violin_cluster',i,'.png'),width = 20,height = 15)
    ggsave(p1,filename = paste0(outdir,'/violin/violin_cluster',i,'.pdf'),width = 20,height = 15)
    p1<-FeaturePlot(integrated_data,features = as.character(top12[top12$cluster==i,]$gene),label = T,reduction = reduction,raster = F) & scale_colour_gradientn(colours = c("lightblue","orange","red"))
    ggsave(p1,filename = paste0(outdir,'/featureplot/featureplot_cluster',i,'.png'),width = 20,height = 15)
    ggsave(p1,filename = paste0(outdir,'/featureplot/featureplot_cluster',i,'.pdf'),width = 20,height = 15)
}
ggsave(p2,filename = paste0(outdir,'/heatmap/top10_heatmap.png'),width = 25,height = 20)
