# .libPaths('/share/apps/R/R-4.1.2-4/lib64/R/library')  # 已注释，使用系统默认R库路径
suppressMessages({
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(Seurat)
  library(argparser)})

argv <- arg_parser('')
argv <- add_argument(argv,"--rds1", help="the merge first seurat rds file")
argv <- add_argument(argv,"--rds2", help="the subcluster seurat rds file")
argv <- add_argument(argv,"--group", help="the group flag")
argv <- add_argument(argv,"--celltype", help="the celltype flag")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- parse_args(argv)

rds1 <- argv$rds1
rds2 <- argv$rds2
group <- argv$group
celltype <- argv$celltype
#subcluster <- argv$subcluster
outdir <- argv$outdir

if(is.na(rds2)){
  each_type_rds<-readRDS(file = rds1)
}
if(!is.na(rds2)){
  original_rds<-readRDS(file = rds1)
  each_type_rds<-readRDS(file = rds2)
}
if(is.na(group)){
  group<-'group'
}
if(is.na(celltype)){
  celltype<-'seurat_clusters'
}
#########################################################
#color_panel<-c('#FFC1C1','#FFA500','#FFF68F','#7FFFD4','#B0E2FF','#FF83FA','#F08080','#54FF9F','#C0FF3E','#6495ED','#8470FF','#FF69B4',
#               '#FF1493','#FFDEAD','#0000FF','#8B1C62','#FF0000','#FFFF00','#00EE00','#B22222','#00BFFF','#006400','#000080','#CD5C5C',
#               '#9400D3','#836FFF','#008B8B','#CD5B45','#2E8B57','#473C8B','#BCD2EE','#363636','#F5F5DC','#53868B','#D8BFD8','#8B5A00',
#               '#8B658B','#C1CDC1','#EED5D2','#000000','#828282','#BC8F8F','#FF4040','#98F5FF','#90EE90')
color_panel <-c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
             '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
             '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
             '#CCC9E6','#3A6963','#68A180','#476D87','#9FA3A8','#968175')
color_panel <-c('#53A85F','#58A4C3','#AB3282','#8C549C','#BD956A','#57C3F3','#6778AE','#F3B1A0','#F1BB72',
             '#DCC1DD','#E95C59','#625D9E','#F7F398','#E63863','#5F3D69','#C5DEBA','#CCE0F5','#B53E2B',
             '#AA9A59','#E39A35','#91D0BE','#23452F','#E4C755','#585658','#C1E6F3','#D6E7A3','#712820',
             '#CCC9E6','#3A6963','#68A180','#61B934','#476D87','#9FA3A8','#968175','#BC8F8F','#D8BFD8')
plot_cell_number_group1 = function(all_cluster_data){
  p1 = ggplot(all_cluster_data,mapping = aes(x=group, y=Proportion, fill=cell_type))+
    geom_bar(stat = 'identity',width = 0.2,position = "fill", col='black')+
    labs(y='proportion',x='group',fill="cluster")+ 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12,color = "black", hjust=1,vjust=1,angle=45),
          axis.line = element_line(colour = "black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=color_panel)#+coord_flip()
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_group.png'),width = 10,height = 8)
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_group.pdf'),width = 10,height = 8)
  
  p2 = ggplot(all_cluster_data, mapping = aes(x=cell_type, y=Proportion, fill=group))+
    facet_wrap(~cell_type,scales = 'free')+
    geom_bar(stat = 'identity',position=position_dodge(0.8),width=0.7)+
    scale_y_continuous(expand = c(0,0))+
    labs(x="",y="Percentage(%)")+ theme_bw()+
    theme(axis.title.y = element_text(size=15,color = "black"),
          axis.text = element_text(size=15,color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 10,color="black"),
          legend.text = element_text(size = 10,color="black"))
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_group_facet.png'),width = 15,height = 12)
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_group_facet.pdf'),width = 15,height = 12)
  
  p3 = ggplot(all_cluster_data, mapping = aes(x=group, y=cell_number, fill=group))+
    facet_wrap(~cell_type,scales = 'free_y')+
    geom_bar(stat = 'identity',position=position_dodge(0.8),width=0.7)+
        geom_text(aes(label = cell_number),
          position = position_dodge(width = 0.7),
          vjust = -0.2,
          size = 4)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    labs(x="",y="Cell number")+ theme_bw()+
    theme(axis.title.y = element_text(size=15,color = "black"),
          axis.text = element_text(size=15,color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 10,color="black"),
          legend.text = element_text(size = 10,color="black"))
  ggsave(p3,filename = paste0(outdir,'/cell_number_group_facet.png'),width = 15,height = 12)
  ggsave(p3,filename = paste0(outdir,'/cell_number_group_facet.pdf'),width = 15,height = 12)

}
plot_cell_number_group2 = function(all_cluster_data){
  p1 = ggplot(all_cluster_data,mapping = aes(x=group, y=Proportion, fill=cell_type))+
    geom_bar(stat = 'identity',width = 0.2, col='black')+
    labs(y='proportion',x='group',fill="cluster")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12,color = "black", hjust=1,vjust=1,angle=45),
          axis.line = element_line(colour = "black"))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=color_panel)#+coord_flip()
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_group.png'),width = 10,height = 8)
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_group.pdf'),width = 10,height = 8)

  p2 = ggplot(all_cluster_data, mapping = aes(x=cell_type, y=Proportion, fill=group))+
    facet_wrap(~cell_type,scales = 'free')+
    geom_bar(stat = 'identity',position=position_dodge(0.8),width=0.7)+
    scale_y_continuous(expand = c(0,0))+
    labs(x="",y="Percentage(%)")+ theme_bw()+
    theme(axis.title.y = element_text(size=15,color = "black"),
          axis.text = element_text(size=15,color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 10,color="black"),
          legend.text = element_text(size = 10,color="black"))
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_group_facet.png'),width = 15,height = 12)
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_group_facet.pdf'),width = 15,height = 12)

  p3 = ggplot(all_cluster_data, mapping = aes(x=group, y=cell_number, fill=group))+
    facet_wrap(~cell_type,scales = 'free_y')+
    geom_bar(stat = 'identity',position=position_dodge(0.8),width=0.7)+
    geom_text(aes(label = cell_number),
          position = position_dodge(width = 0.7),
          vjust = -0.2,
          size = 4)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    labs(x="",y="Cell number")+ theme_bw()+
    theme(axis.title.y = element_text(size=15,color = "black"),
          axis.text = element_text(size=15,color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 10,color="black"),
          legend.text = element_text(size = 10,color="black"))
  ggsave(p3,filename = paste0(outdir,'/cell_number_group_facet.png'),width = 15,height = 12)
  ggsave(p3,filename = paste0(outdir,'/cell_number_group_facet.pdf'),width = 15,height = 12)

}
plot_cell_number_sample1 = function(all_cluster_data){
  p1<-ggplot(all_cluster_data,mapping = aes(x=sample, y=Proportion,fill=cell_type))+
    geom_bar(stat = 'identity',width = 0.2,position = "fill", col='black')+
    labs(y='proportion',x='sample',fill="cluster")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12,color = "black", hjust=1,vjust=1,angle=45),
          axis.line = element_line(colour = "black"))+
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=color_panel)#+ coord_flip()
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_sample.png'),width = 12,height = 8)
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_sample.pdf'),width = 12,height = 8)
  
  p2 = ggplot(all_cluster_data, mapping = aes(x=cell_type, y=Proportion, fill=sample))+
    facet_wrap(~cell_type,scales = 'free')+
    geom_bar(stat = 'identity',position=position_dodge(0.8),width=0.7)+
    scale_y_continuous(expand = c(0,0))+
    labs(x="",y="Percentage(%)")+ theme_bw()+
    theme(axis.title.y = element_text(size=15,color = "black"),
          axis.text = element_text(size=15,color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 10,color="black"),
          legend.text = element_text(size = 10,color="black"))
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_sample_facet.png'),width = 18,height = 15)
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_sample_facet.pdf'),width = 18,height = 15)
}
plot_cell_number_sample2 = function(all_cluster_data){
  p1<-ggplot(all_cluster_data,mapping = aes(x=sample, y=Proportion,fill=cell_type))+
    geom_bar(stat = 'identity',width = 0.2, col='black')+
    labs(y='proportion',x='sample',fill="cluster")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12,color = "black", hjust=1,vjust=1,angle=45),
          axis.line = element_line(colour = "black"))+
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=color_panel)#+ coord_flip()
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_sample.png'),width = 12,height = 8)
  ggsave(p1,filename = paste0(outdir,'/cell_number_rate_sample.pdf'),width = 12,height = 8)

  p2 = ggplot(all_cluster_data, mapping = aes(x=cell_type, y=Proportion, fill=sample))+
    facet_wrap(~cell_type,scales = 'free')+
    geom_bar(stat = 'identity',position=position_dodge(0.8),width=0.7)+
    scale_y_continuous(expand = c(0,0))+
    labs(x="",y="Percentage(%)")+ theme_bw()+
    theme(axis.title.y = element_text(size=15,color = "black"),
          axis.text = element_text(size=15,color = "black"),
          axis.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15,color="black"),
          legend.title = element_text(size = 10,color="black"),
          legend.text = element_text(size = 10,color="black"))
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_sample_facet.png'),width = 18,height = 15)
  ggsave(p2,filename = paste0(outdir,'/cell_number_rate_sample_facet.pdf'),width = 18,height = 15)
}
#########################################################
####### cellnumber_stat
print("create the cell_number_stat file!")
all_df = data.frame()
each_type_rds$group = factor(each_type_rds@meta.data[,group])
each_type_rds$celltype<-factor(each_type_rds@meta.data[,celltype])
each_type_rds$orig.ident = factor(each_type_rds$orig.ident)
all_df = data.frame(table(each_type_rds@meta.data$orig.ident))
colnames(all_df) = c("sample","cell_number")
write.table(all_df,file = paste0(outdir, '/cell_number_stat.xls'),quote = F,sep = '\t',row.names = F)

####### without subcluster
print("Cell Ratio Section of Cluster!")
if(is.na(rds2)){
  ####### plot by group
  all_cluster_data = data.frame(); 
  for (i in as.vector(levels(each_type_rds@meta.data$celltype))) {
    for (each in as.vector(levels(each_type_rds@meta.data$group))) {
      each_cluster_data<-data.frame(cluster=i,group=each,cell_number=table(each_type_rds@meta.data[each_type_rds@meta.data$group==each,]$celltype)[[i]],
                                    rate=round(100*table(each_type_rds@meta.data[each_type_rds@meta.data$group==each,]$celltype)[[i]]/table(each_type_rds@meta.data$group)[[each]],digits = 2))
      all_cluster_data<-rbind(all_cluster_data,each_cluster_data)
    }
  }
  write.table(all_cluster_data,file = paste0(outdir,'/cell_number_group.xls'),quote = F,sep = '\t',row.names = F)
  colnames(all_cluster_data)<-c('cell_type','group','cell_number','Proportion')
  all_cluster_data$cell_type<-factor(all_cluster_data$cell_type,levels=levels(each_type_rds@meta.data$celltype))
  all_cluster_data$group<-factor(all_cluster_data$group,levels=levels(each_type_rds@meta.data$group))
  plot_cell_number_group1(all_cluster_data)
  
  ####### plot by sample
  all_cluster_data = data.frame(); 
  for (i in as.vector(levels(each_type_rds@meta.data$celltype))) {
    for (each in as.vector(levels(each_type_rds@meta.data$orig.ident))) {
      each_cluster_data<-data.frame(cluster=i,sample=each,
                                    cell_number=table(each_type_rds@meta.data[each_type_rds@meta.data$orig.ident==each,]$celltype)[[i]],
                                    rate=round(100*table(each_type_rds@meta.data[each_type_rds@meta.data$orig.ident==each,]$celltype)[[i]]/table(each_type_rds@meta.data$orig.ident)[[each]],digits = 2))
      all_cluster_data<-rbind(all_cluster_data,each_cluster_data)
    }
  }
  write.table(all_cluster_data,file = paste0(outdir,'/cell_number_sample.xls'),quote = F,sep = '\t',row.names = F)
  colnames(all_cluster_data)<-c('cell_type','sample','cell_number','Proportion')
  all_cluster_data$cell_type<-factor(all_cluster_data$cell_type,levels=levels(each_type_rds@meta.data$celltype))
  all_cluster_data$sample<-factor(all_cluster_data$sample,levels=levels(each_type_rds@meta.data$orig.ident))
  plot_cell_number_sample1(all_cluster_data)

  ####### plot by cluster
  all_cluster_data = data.frame();
  for (i in as.vector(levels(each_type_rds@meta.data$celltype))) {
    for (each in as.vector(levels(each_type_rds@meta.data$group))) {
      each_cluster_data<-data.frame(cluster=i,group=each,cell_number=table(each_type_rds@meta.data[each_type_rds@meta.data$celltype==i,]$group)[[each]],
                                  rate=round(100*table(each_type_rds@meta.data[each_type_rds@meta.data$celltype==i,]$group)[[each]]/table(each_type_rds@meta.data$celltype)[[i]],digits = 2))
    all_cluster_data<-rbind(all_cluster_data,each_cluster_data)
    }
  }
  write.table(all_cluster_data,file = paste0(outdir,'/cell_number_cluster.xls'),quote = F,sep = '\t',row.names = F)
  colnames(all_cluster_data)<-c('cell_type','group','cell_number','Proportion')
  all_cluster_data$cell_type<-factor(all_cluster_data$cell_type,levels=levels(each_type_rds@meta.data$celltype))
  all_cluster_data$group<-factor(all_cluster_data$group,levels=levels(each_type_rds@meta.data$group))
  p<-ggplot(all_cluster_data,mapping = aes(x=cell_type, y=Proportion,fill=group))+geom_bar(stat = 'identity',width = 0.2)+labs(y='proportion(%)',x='cluster',fill="group")+ theme(panel.grid.major = element_blank(),axis.text.x = element_text(size=12,color = "black", hjust=1,vjust=1,angle=45), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  ggsave(p,filename = paste0(outdir,'/cell_number_rate_cluster.png'),width = 12,height = 8)
  ggsave(p,filename = paste0(outdir,'/cell_number_rate_cluster.pdf'),width = 12,height = 8)

}

####### with subcluster
print("Cell Ratio Section of Subcluster!")
if(!is.na(rds2)){
  original_rds$group = factor(original_rds@meta.data[,group])
  original_rds$celltype<-factor(original_rds@meta.data[,celltype])
  ####### plot by group
  all_cluster_data = data.frame(); 
  for (i in as.vector(levels(each_type_rds@meta.data$celltype))) {
    for (each in as.vector(levels(each_type_rds@meta.data$group))) {
      each_cluster_data<-data.frame(cluster=i,group=each,cell_number=table(each_type_rds@meta.data[each_type_rds@meta.data$group==each,]$celltype)[[i]],
                                    rate=round(100*table(each_type_rds@meta.data[each_type_rds@meta.data$group==each,]$celltype)[[i]]/table(original_rds@meta.data$group)[[each]],digits = 2))
      all_cluster_data<-rbind(all_cluster_data,each_cluster_data)
    }
  }
  write.table(all_cluster_data,file = paste0(outdir,'/cell_number_group.xls'),quote = F,sep = '\t',row.names = F)
  colnames(all_cluster_data)<-c('cell_type','group','cell_number','Proportion')
  all_cluster_data$cell_type<-factor(all_cluster_data$cell_type,levels=levels(each_type_rds@meta.data$celltype))
  all_cluster_data$group<-factor(all_cluster_data$group,levels=levels(each_type_rds@meta.data$group))
  plot_cell_number_group2(all_cluster_data)
  
  ####### plot by sample
  all_cluster_data = data.frame(); 
  for (i in as.vector(levels(each_type_rds@meta.data$celltype))) {
    for (each in as.vector(levels(each_type_rds@meta.data$orig.ident))) {
      each_cluster_data<-data.frame(cluster=i,sample=each,
                                    cell_number=table(each_type_rds@meta.data[each_type_rds@meta.data$orig.ident==each,]$celltype)[[i]],
                                    rate=round(100*table(each_type_rds@meta.data[each_type_rds@meta.data$orig.ident==each,]$celltype)[[i]]/table(original_rds@meta.data$orig.ident)[[each]],digits = 2))
      all_cluster_data<-rbind(all_cluster_data,each_cluster_data)
    }
  }
  write.table(all_cluster_data,file = paste0(outdir,'/cell_number_sample.xls'),quote = F,sep = '\t',row.names = F)
  colnames(all_cluster_data)<-c('cell_type','sample','cell_number','Proportion')
  all_cluster_data$cell_type<-factor(all_cluster_data$cell_type,levels=levels(each_type_rds@meta.data$celltype))
  all_cluster_data$sample<-factor( all_cluster_data$sample,levels=levels(each_type_rds@meta.data$orig.ident))
  plot_cell_number_sample2(all_cluster_data)

  ####### plot by cluster
  all_cluster_data = data.frame();
  for (each in as.vector(levels(each_type_rds@meta.data$group))) {
    for (i in as.vector(levels(each_type_rds@meta.data$celltype))) {
      each_cluster_data<-data.frame(cluster=i,group=each,cell_number=table(each_type_rds@meta.data[each_type_rds@meta.data$celltype==i,]$group)[[each]],rate=round(100*table(each_type_rds@meta.data[each_type_rds@meta.data$celltype==i,]$group)[[each]]/table(each_type_rds@meta.data$celltype)[[i]],digits = 2))
    all_cluster_data<-rbind(all_cluster_data,each_cluster_data)
    }
  }
  write.table(all_cluster_data,file = paste0(outdir,'/cell_number_cluster.xls'),quote = F,sep = '\t',row.names = F)
  colnames(all_cluster_data)<-c('cell_type','group','cell_number','Proportion')
  all_cluster_data$cell_type<-factor(all_cluster_data$cell_type,levels=levels(each_type_rds@meta.data$celltype))
  all_cluster_data$group<-factor(all_cluster_data$group,levels=levels(each_type_rds@meta.data$group))
  p<-ggplot(all_cluster_data,mapping = aes(x=cell_type, y=Proportion,fill=group))+geom_bar(stat = 'identity',width = 0.2)+labs(y='proportion(%)',x='cluster',fill="group")+ theme(panel.grid.major = element_blank(), axis.text.x = element_text(size=12,color = "black", hjust=1,vjust=1,angle=45),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  ggsave(p,filename = paste0(outdir,'/cell_number_rate_cluster.png'),width = 12,height = 8)
  ggsave(p,filename = paste0(outdir,'/cell_number_rate_cluster.pdf'),width = 12,height = 8)

}



