# .libPaths(c("/share/apps/R/R-4.1.2-4/lib64/R/library/"))  # 已注释，使用系统默认R库路径
suppressMessages({
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(argparser)
library(stringr)})
library(clusterProfiler)
library(org.Mm.eg.db)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="the seurat rds file")
argv <- add_argument(argv,"--group", help="the compare of diff,such as AvsB1:AvsB2,AvsB")
argv <- add_argument(argv,"--group_flag", help="the group flag of seurat,such as group,group2")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--species", help="the species of project")
argv <- add_argument(argv,"--volcano_flag", help="the choose volcane genelist flag,1 or 0, 1:all gene,0:partgene")
argv <- add_argument(argv,"--p_val_adj", help="the value of p_val_adj in volcane")
argv <- add_argument(argv,"--avg_log2FC", help="the value of avg_log2FC in volcane")
argv <- parse_args(argv)

library(ggrepel)
Run_diffgene_volcano<-function(diffgene, label_gene_list, cpname, p_val_adj, logFC, outdir){
    volcano <- subset(diffgene,select=c(gene,avg_log2FC,p_val_adj))
    min_p_val_adj<-ceiling(-log10(min(volcano[volcano$p_val_adj!=0,]$p_val_adj)))
    volcano$log_p_val_adj<--log10(volcano$p_val_adj)
    if('Inf' %in% volcano$log_p_val_adj){
      volcano[volcano$log_p_val_adj=='Inf',]$log_p_val_adj<-min_p_val_adj+100}
    p_val_adj<-as.numeric(p_val_adj)
    logFC<-as.numeric(logFC)
    mark_label<-c()
    for(each in volcano$gene){
      if(each %in% label_gene_list){
       mark_label<-c(mark_label,each)
      }else{
      mark_label<-c(mark_label,"")
      }
    }
    volcano$label<-mark_label
    volcano["group"] <- "NO"
    avg_log2FC=logFC
    avg_log2FC2=-logFC
    volcano[which(volcano["p_val_adj"] <= p_val_adj & volcano["avg_log2FC"] >= avg_log2FC),"group"] <- "UP"
    volcano[which(volcano["p_val_adj"] <= p_val_adj & volcano["avg_log2FC"] <= avg_log2FC2),"group"] <- "DOWN"
    NO_number <- nrow(volcano[which(volcano["group"]=='NO'),])
    UP_number <- nrow(volcano[which(volcano["group"]=='UP'),])
    DOWN_number <- nrow(volcano[which(volcano["group"]=='DOWN'),])
    volcano[which(volcano["group"]=='NO'),"group"] <- paste('NO',NO_number,sep=' ')
    volcano[which(volcano["group"]=='UP'),"group"] <- paste('UP',UP_number,sep=' ')
    volcano[which(volcano["group"]=='DOWN'),"group"] <- paste('DOWN',DOWN_number,sep=' ')
    volcano$group <- factor(volcano$group,levels=c(paste('UP',UP_number,sep=' '),paste('DOWN',DOWN_number,sep=' '),paste('NO',NO_number,sep=' ')))
    fc_min <- min(volcano$avg_log2FC)
    fc_max <- max(volcano$avg_log2FC)
    x_breaks <- c(-log10(p_val_adj),pretty(-log10(volcano$p_val_adj),n=3))
    if (length(unique(round(pretty(-log10(volcano$p_val_adj),n=3))))==length(round(pretty(-log10(volcano$p_val_adj),n=3)))){
    x_labels <- c(round(-log10(p_val_adj),3),round(pretty(-log10(volcano$p_val_adj),n=3)))}else{
    x_labels <- c(round(-log10(p_val_adj),3),pretty(-log10(volcano$p_val_adj),n=3))
    }
    P <- ggplot(volcano,aes(x=avg_log2FC,y=log_p_val_adj,group=group,color=group)) +
       geom_point(size=4) +ylab('-log10(p_val_adj)')+
       labs(title=cpname,color=paste('p_val_adj<',p_val_adj,sep='')) +
       #geom_text(aes(label = label), size = 3,vjust=-0.5, alpha=0.8)+
       theme(plot.title=element_text(hjust=0.5),text = element_text(size = 20)) +
       scale_x_continuous(breaks=c(-logFC,logFC,round(seq(fc_min,fc_max,by=(fc_max-fc_min)/4)))) +
       scale_y_continuous(breaks=x_breaks[x_breaks!=0],labels=x_labels[x_labels!=0]) +
       geom_vline(xintercept=c(-logFC,logFC),linetype='dotdash',size=0.8,color='grey') +
       geom_hline(yintercept=-log10(p_val_adj),linetype='dotdash',size=0.8,color='grey')+
       geom_text_repel(data = volcano, aes(x = avg_log2FC, y = log_p_val_adj,label = label),size = 3,max.overlaps=160,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))
       #geom_text_repel(data = volcano, aes(x = avg_log2FC, y = log_p_val_adj,label = label),size = 4,box.padding = unit(0.5, "lines"),point.padding = unit(0, "lines"),segment.color = "black",show.legend = F,nudge_y = 0)
       #geom_text_repel(data = volcano, aes(x = avg_log2FC, y = log_p_val_adj,label = label),size = 2.5,box.padding = unit(0.5, "lines"),point.padding = unit(0.5, "lines"),segment.color = "black",show.legend = F,nudge_y = 0.1)
    P <- P + theme(panel.background=element_rect(fill="transparent"),axis.line=element_line())
    plot(P)
    ggsave(P,file=paste(outdir,'/',cpname,'_volcano.png',sep=''),width = 10,height = 10)
    ggsave(P,file=paste(outdir,'/',cpname,'_volcano.pdf',sep=''),width = 10,height = 10)
    return(volcano)
}
go_enrich_analysis_all<-function(genes,species,outdir,prefix){
  if(species=='Human'){
    #library(org.Hs.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Hs.eg.db'
  }else if(species=='Mouse'){
    #library(org.Mm.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Mm.eg.db'
  }else if(species=='Rat'){
    #library(org.Rn.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Rn.eg.db'
  }
  GainALLGo<-enrichGO(genes, db, keyType ='SYMBOL', ont="ALL", pvalueCutoff=1,qvalueCutoff = 1)
  if(class(GainALLGo)!="NULL"){
  GOenrich <- as.data.frame(GainALLGo)
  #GOenrich <- data.frame(GOenrich[,1:6],GOenrich[8],GOenrich[9])
  GOenrich <- data.frame(GOenrich)
  write.table(GOenrich, file=paste0(outdir,"/",prefix,"_GO.ALL.txt"), sep="\t", quote=F, row.names=F)
  all_bar<-barplot(GainALLGo, showCategory=15, split="ONTOLOGY",label_format=100)+ 
               theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))+
          facet_grid(ONTOLOGY~., scale='free')
  ggsave(all_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentALL.bar.png"),dpi=400, width = 18, height = 20)
  ggsave(all_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentALL.bar.pdf"),width = 18, height = 20)
}}

go_enrich_analysis_all_dot<-function(genes,species,outdir,prefix){
  if(species=='Human'){
    #library(org.Hs.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Hs.eg.db'
  }else if(species=='Mouse'){
    #library(org.Mm.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Mm.eg.db'
  }else if(species=='Rat'){
    #library(org.Rn.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Rn.eg.db'
  }
  GainALLGo<-enrichGO(genes, db, keyType ='SYMBOL', ont="ALL", pvalueCutoff=1,qvalueCutoff = 1)
  if(class(GainALLGo)!="NULL"){
  GOenrich <- as.data.frame(GainALLGo)
  #GOenrich <- data.frame(GOenrich[,1:6],GOenrich[8],GOenrich[9])
  GOenrich <- data.frame(GOenrich)
  write.table(GOenrich, file=paste0(outdir,"/",prefix,"_GO.ALL.txt"), sep="\t", quote=F, row.names=F)
  all_bar<-dotplot(GainALLGo, showCategory=15, split="ONTOLOGY",label_format=100)+
               theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))+
          facet_grid(ONTOLOGY~., scale='free')
  ggsave(all_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentALL.dot.png"),dpi=400, width = 18, height = 20)
  ggsave(all_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentALL.dot.pdf"),width = 18, height = 20)
}}

go_enrich_analysis<-function(genes,species,outdir,prefix){
  if(species=='Human'){
    #library(org.Hs.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Hs.eg.db'
  }else if(species=='Mouse'){
    #library(org.Mm.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Mm.eg.db'
  }else if(species=='Rat'){
    #library(org.Rn.eg.db,lib.loc = '/share/apps/R/R-4.2.0/lib64/R/library')
    db='org.Rn.eg.db'
  }
  #BP
  GainBPGo<-enrichGO(genes, db, keyType ='SYMBOL', ont="BP", pvalueCutoff=1,qvalueCutoff = 1)
  if(class(GainBPGo)!="NULL"){
    GOenrich <- as.data.frame(GainBPGo)
    GOenrich <- data.frame(GOenrich[,1:6],GOenrich[8],GOenrich[9])
    write.table(GOenrich, file=paste0(outdir,"/",prefix,"_GO.BP.txt"), sep="\t", quote=F, row.names=F)
    bp_bar<-barplot(GainBPGo, showCategory=30,label_format=100)+
    theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))
    ggsave(bp_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentBP.bar.png"),dpi=400, width = 15, height = 15)
    ggsave(bp_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentBP.bar.pdf"),width = 15, height = 15)
    #pdf(paste(outdir,"/",prefix,".GObp_DAG.pdf",sep=""))
    #plotGOgraph(GainBPGo,firstSigNodes = 5)
    #dev.off()
    #png(paste(outdir,"/",prefix,".GObp_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4)
    #plotGOgraph(GainBPGo,firstSigNodes = 5)
    #dev.off()
  }
  #MF
  GainMFGo<-enrichGO(genes, db, keyType ='SYMBOL', ont="MF", pvalueCutoff=1,qvalueCutoff = 1)
  if(class(GainMFGo)!="NULL"){
    GOenrich <- as.data.frame(GainMFGo)
    GOenrich <- data.frame(GOenrich[,1:6],GOenrich[8],GOenrich[9])
    write.table(GOenrich, file=paste0(outdir,"/",prefix,"_GO.MF.txt"), sep="\t", quote=F, row.names=F)
    mf_bar<-barplot(GainMFGo, showCategory=30,label_format=100)+
    theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))
    ggsave(mf_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentMF.bar.png"),dpi=400, width = 15, height = 15)
    ggsave(mf_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentMF.bar.pdf"),width = 15, height = 15)
    #pdf(paste(outdir,"/",prefix,".GOmf_DAG.pdf",sep=""))
    #plotGOgraph(GainMFGo,firstSigNodes = 5)
    #dev.off()
    #png(paste(outdir,"/",prefix,".GOmf_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4)
    #plotGOgraph(GainMFGo,firstSigNodes = 5)
    #dev.off()
  }
  #CC
  GainCCGo<-enrichGO(genes, db, keyType ='SYMBOL', ont="CC", pvalueCutoff=1,qvalueCutoff = 1)
  if(class(GainCCGo)!="NULL"){
    GOenrich <- as.data.frame(GainCCGo)
    GOenrich <- data.frame(GOenrich[,1:6],GOenrich[8],GOenrich[9])
    write.table(GOenrich, file=paste0(outdir,"/",prefix,"_GO.CC.txt"), sep="\t", quote=F, row.names=F)
    cc_bar<-barplot(GainCCGo, showCategory=30,label_format=100)+
    theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))
    ggsave(cc_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentCC.bar.png"),dpi=400, width = 15, height = 15)
    ggsave(cc_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentCC.bar.pdf"),width = 15, height = 15)
    #pdf(paste(outdir,"/",prefix,".GOcc_DAG.pdf",sep=""))
    #plotGOgraph(GainCCGo,firstSigNodes = 5)
    #dev.off()
    #png(paste(outdir,"/",prefix,".GOcc_DAG.png",sep=""),type="cairo-png",width=480*4,height=480*4,res=72*4)
    #plotGOgraph(GainCCGo,firstSigNodes = 5)
    #dev.off()
  }
}

kegg_enrich_analysis<-function(genes,species,outdir,prefix){
  if(species=='Human'){
    entrezIDs <- unlist(mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA))
    enrich_kegg<-enrichKEGG(gene = entrezIDs, organism = "hsa", keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", universe = keys(org.Hs.eg.db), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = T)
  }else if(species=='Mouse'){
    entrezIDs <- unlist(mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA))
    enrich_kegg<-enrichKEGG(gene = entrezIDs, organism = "mmu", keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", universe = keys(org.Mm.eg.db), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = T)
  }else if(species=='Rat'){
    entrezIDs <- unlist(mget(genes, org.Rn.egSYMBOL2EG, ifnotfound=NA))
    enrich_kegg<-enrichKEGG(gene = entrezIDs, organism = "rno", keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", universe = keys(org.Rn.eg.db), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = T)
  }
  if(class(enrich_kegg)!="NULL"){
    gene_list <- strsplit(enrich_kegg$geneID,split='/')
    geneName <- unlist(lapply(gene_list,FUN=function(x){paste(names(entrezIDs[match(x,entrezIDs)]),collapse='/')}))
    KEGGenrich <- as.data.frame(enrich_kegg)
    KEGGenrich <- data.frame(KEGGenrich[,1:6],KEGGenrich[8],geneName,KEGGenrich[9])
    write.table(KEGGenrich ,file=paste0(outdir,"/",prefix,"_KEGG.txt"),sep="\t",quote=F,row.names = F)
    kegg_bar<-barplot(enrich_kegg, showCategory=30,label_format=100)+
    theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))
    ggsave(kegg_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentKEGG.bar.png"),dpi=400, width = 15, height = 15)
    ggsave(kegg_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentKEGG.bar.pdf"),width = 15, height = 15)}
}

kegg_enrich_analysis_dot<-function(genes,species,outdir,prefix){
  if(species=='Human'){
    entrezIDs <- unlist(mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA))
    enrich_kegg<-enrichKEGG(gene = entrezIDs, organism = "hsa", keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", universe = keys(org.Hs.eg.db), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = T)
  }else if(species=='Mouse'){
    entrezIDs <- unlist(mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA))
    enrich_kegg<-enrichKEGG(gene = entrezIDs, organism = "mmu", keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", universe = keys(org.Mm.eg.db), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = T)
  }else if(species=='Rat'){
    entrezIDs <- unlist(mget(genes, org.Rn.egSYMBOL2EG, ifnotfound=NA))
    enrich_kegg<-enrichKEGG(gene = entrezIDs, organism = "rno", keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", universe = keys(org.Rn.eg.db), minGSSize = 10, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = T)
  }
  if(class(enrich_kegg)!="NULL"){
    gene_list <- strsplit(enrich_kegg$geneID,split='/')
    geneName <- unlist(lapply(gene_list,FUN=function(x){paste(names(entrezIDs[match(x,entrezIDs)]),collapse='/')}))
    KEGGenrich <- as.data.frame(enrich_kegg)
    KEGGenrich <- data.frame(KEGGenrich[,1:6],KEGGenrich[8],geneName,KEGGenrich[9])
    write.table(KEGGenrich ,file=paste0(outdir,"/",prefix,"_KEGG.txt"),sep="\t",quote=F,row.names = F)
    kegg_bar<-dotplot(enrich_kegg, showCategory=30,label_format=100)+
    theme(axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14,colour="black"),
          axis.title.x.bottom = element_text(size=14,colour="black"))
    ggsave(kegg_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentKEGG.dot.png"),dpi=400, width = 15, height = 15)
    ggsave(kegg_bar, filename = paste0(outdir,"/",prefix,"_EnrichmentKEGG.dot.pdf"),width = 15, height = 15)}
}

rds <- argv$rds
group <- argv$group
group_flag <- argv$group_flag
outdir <- argv$outdir
species <- argv$species
volcano_flag <- argv$volcano_flag
p_val_adj <- argv$p_val_adj
avg_log2FC <- argv$avg_log2FC

if(is.na(p_val_adj)){
  p_val_adj = 0.05
}else{p_val_adj = as.numeric(p_val_adj)}
if(is.na(avg_log2FC)){
  avg_log2FC = 0
}else{avg_log2FC = as.numeric(avg_log2FC)}

group<-str_split(group,pattern = ',')[[1]]
group_flag<-str_split(group_flag,pattern = ',')[[1]]
pancreas.integrated<-readRDS(rds)

DefaultAssay(pancreas.integrated)<-'RNA'
Idents(pancreas.integrated)<-pancreas.integrated@meta.data$seurat_clusters

pos = 1
for(each_flag in group_flag){
  each_group<-str_split(group[pos],pattern = ':')[[1]]
  for(each in each_group){
    dir.create(paste0(outdir,'/diff/'))
    dir.create(paste0(outdir,'/diff/',each))
    dir.create(paste0(outdir,'/diff/',each,'/volcano'))
    dir.create(paste0(outdir,'/enrich/'))
    dir.create(paste0(outdir,'/enrich/',each))
    dir.create(paste0(outdir,'/enrich/',each,'/ALL'))
    dir.create(paste0(outdir,'/enrich/',each,'/UP'))
    dir.create(paste0(outdir,'/enrich/',each,'/DOWN'))
    diff1 <- str_split(each,pattern = 'vs')[[1]][1]
    diff2 <- str_split(each,pattern = 'vs')[[1]][2]
    for (i in as.vector(levels(pancreas.integrated@meta.data$seurat_clusters))) {
      if(table(pancreas.integrated@meta.data[pancreas.integrated@meta.data$seurat_clusters==i,][,each_flag])[[diff1]]>=3 & table(pancreas.integrated@meta.data[pancreas.integrated@meta.data$seurat_clusters==i,][,each_flag])[[diff2]]>=3){
      dir.create(paste0(outdir,'/enrich/',each,'/ALL/cluster',i))
      dir.create(paste0(outdir,'/enrich/',each,'/UP/cluster',i))
      dir.create(paste0(outdir,'/enrich/',each,'/DOWN/cluster',i))
pancreas.integrated <- JoinLayers(pancreas.integrated)
      pancreas.integrated.marker <- FindMarkers(pancreas.integrated, ident.1 = diff1,ident.2 = diff2,group.by = each_flag, subset.ident = i,only.pos = FALSE,min.cells.group = 0)
      pancreas.integrated.marker$gene<-rownames(pancreas.integrated.marker)
      pancreas.integrated.marker<-pancreas.integrated.marker[,c(ncol(pancreas.integrated.marker),1:ncol(pancreas.integrated.marker)-1)]
      write.table(pancreas.integrated.marker,paste0(outdir,'/diff/',each,'/cluster',i,'_',each,'_diffgene.xls'),quote = F,row.names = F,sep = '\t')
      sig_mark<-subset(pancreas.integrated.marker,p_val_adj<0.05)
      write.table(sig_mark,paste0(outdir,'/diff/',each,'/cluster',i,'_',each,'_sig_diffgene.xls'),quote = F,row.names = F,sep = '\t')
      if(volcano_flag==1){
        gene_label <- as.character(sig_mark$gene)
        Run_diffgene_volcano(diffgene=pancreas.integrated.marker,label_gene_list=gene_label,cpname=paste0('cluster',i,'_',each),p_val_adj=p_val_adj,logFC=avg_log2FC,outdir=paste0(outdir,'/diff/',each,'/volcano'))}
      if(nrow(sig_mark)>=50 & volcano_flag==0){
        sig_mark$avg_log2FC2<-abs(sig_mark$avg_log2FC)
        top50 <- sig_mark %>% top_n(n = 50, wt = avg_log2FC2)
        gene_label <- as.character(top50$gene)
        Run_diffgene_volcano(diffgene=pancreas.integrated.marker,label_gene_list=gene_label,cpname=paste0('cluster',i,'_',each),p_val_adj=p_val_adj,logFC=avg_log2FC,outdir=paste0(outdir,'/diff/',each,'/volcano'))}
      if(nrow(sig_mark)<50 & volcano_flag==0){
        gene_label <- as.character(sig_mark$gene)
        Run_diffgene_volcano(diffgene=pancreas.integrated.marker,label_gene_list=gene_label,cpname=paste0('cluster',i,'_',each),p_val_adj=p_val_adj,logFC=avg_log2FC,outdir=paste0(outdir,'/diff/',each,'/volcano'))}
      if(!(is.na(species))){
        sig_gene<-sig_mark$gene
        outdir2=paste0(outdir,'/enrich/',each,'/ALL/cluster',i)
        prefix=paste0('cluster',i)
        go_enrich_analysis(sig_gene,species,outdir2,prefix)
        kegg_enrich_analysis(sig_gene,species,outdir2,prefix)
        sig_mark_up<-subset(sig_mark,avg_log2FC>0)
        sig_mark_down<-subset(sig_mark,avg_log2FC<0)
        sig_gene_up<-sig_mark_up$gene
        outdir2=paste0(outdir,'/enrich/',each,'/UP/cluster',i)
        prefix=paste0('cluster',i,'_up')
        go_enrich_analysis(sig_gene_up,species,outdir2,prefix)
        kegg_enrich_analysis(sig_gene_up,species,outdir2,prefix)
        sig_gene_down<-sig_mark_down$gene
        outdir2=paste0(outdir,'/enrich/',each,'/DOWN/cluster',i)
        prefix=paste0('cluster',i,'_down')
        go_enrich_analysis(sig_gene_down,species,outdir2,prefix)
        kegg_enrich_analysis(sig_gene_down,species,outdir2,prefix)
      }}
    }
  }
  pos <-pos+1
}
