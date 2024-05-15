library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
require("RColorBrewer")
library(dplyr)
library(harmony)
library(monocle)

# args<-commandArgs(TRUE)
# cellID<-args[1]

setwd('/Lustre02/zhangJM/01.Project/14.Chicken.Broiler.vs.Layer.study/snRNA.analysis/01.snRNAseq.analysis.for.150g.data.final/Suerat.scale.harmony/Rename.analysis/02.Each.cell.type.for.subCluster/MuSC/Rename.data/Monocle2.cell.trajectory.analysis.test/')


pbmcnewer<-readRDS('merge_processed.anchor.for.rmbatch.rename.rds')


pbmc<-pbmcnewer

data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = pbmc@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

bm280k<-pbmcnewer
# bm280k <- NormalizeData(bm280k, normalization.method = "LogNormalize", scale.factor = 10000)
# bm280k <- FindVariableFeatures(bm280k, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(bm280k)
# bm280k <- ScaleData(bm280k, features = all.genes)
# bm280k <- RunPCA(bm280k,npcs = 50, verbose = FALSE)

# bm280k.integrated<-RunHarmony(bm280k,group.by.vars = 'orig.ident',max.iter.harmony=30,lambda=2)
# bm280k.integrated <- RunUMAP(bm280k.integrated, reduction = "harmony", dims = 1:30)

# p1<-DimPlot(bm280k, reduction = "umap", group.by="integrated_merge_cluster",pt.size = 0.6,label = TRUE, repel = TRUE)
# oup<-paste('UMAP.for.integrated_merge_cluster.new.for.monocle.pdf',sep='')
# ggsave(oup, plot=p1, width = 8, height =8)


deg.cluster <- FindAllMarkers(bm280k,assay='RNA',group.by='integrated_merge_cluster')
write.table(deg.cluster,file=paste('all.cluster.DEgenes.new.txt',sep=''),sep='\t',quote=F,row.names=T,col.names=T)

diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, diff.genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

cds <- orderCells(cds)

p1<-plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(paste('monocle.Pseudotime.UMAP.pdf',sep=""),limitsize=F, plot=p1, width = 8, height = 8)

p2<-plot_cell_trajectory(cds, color_by = "seurat_clusters")
ggsave(paste('monocle.seurat_clusters.UMAP.pdf',sep=""),limitsize=F, plot=p1, width = 8, height = 8)

p3<-plot_cell_trajectory(cds, color_by = "State")
ggsave(paste('monocle.State.UMAP.pdf',sep=""),limitsize=F, plot=p1, width = 8, height = 8)


saveRDS(cds, file = "subcluster.anchor.for.monocle.new.rds")
