#volcano_plot
setwd("~/yinxiaolan")
res <- read.csv("tmp01", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-12,30),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

library(stringr)
library(sva)
library(data.table)
library(readxl)
library(DESeq2)
library(DESeq)
library(pamr)
library(ggpubr)
library(Seurat)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(ggplot2)
library(gplots)
library(pca3d)
library(rgl)
library(scatterplot3d)
library(FactoMineR)
library(ggfortify)
library(useful)
library(tidyverse)
library(kableExtra)
library(xfun)
library(psych)
library(limma)
library(calibrate)
library(pheatmap)
library(apeglm)


#change ensg to geneid and averageif
setwd("~/yinxiaolan")
DF<-read.table("tmp02",header = TRUE)
b<-aggregate(DF[, -c(1:2)], by=list(DF$EntrzID, DF$Name), mean)
write.table(b, file="yinxiaolan_averageif.txt")

#pca 
setwd("~/yinxiaolan")
gse184053<-read.table("tmp03",header=TRUE)
write.csv(gse184053,"gse184053.csv")
colnames(gse184053)
#can not remove batch effect, the results of combat_gse184053 are all NaN.
#bat <-as.data.frame(cbind(colnames(gse184053),c(c(1:9)),c(rep(1, 3),rep(2, 3),rep(3, 3))))
#mod<-model.matrix(~as.factor(bat[,3]), data=bat)
#combat_gse184053<-ComBat(dat = as.matrix(gse184053,),batch = (bat[,2]), par.prior = TRUE)#PCA post combat

#pca
gse184053_01<-gse184053[,-1]
pca<-prcomp(t(gse184053_01))

plot(pca)
{plot(pca$x[,1],pca$x[,2])
  points(pca$x[1:7,1],pca$x[1:7,2], main="Top 2 PCs",col=2)
  points(pca$x[8:14,1],pca$x[8:14,2], main="Top 2 PCs",col=3)
}
which(pca$x[,1]< -50)
which(pca$x[,2]< -50)
percentVar <- pca$sdev^2 / sum( pca$sdev^2)

#3d pca  1. 因为不管怎么筛选，我们都只用到了3个基因的表达量。no prcomp
setwd("~/yinxiaolan")
data4 <- read.table("tmp03", header=T, row.names=1,sep="\t")
data4_nonzero <- data4[rowSums(data4)!=0,]#去除表达值全为0的行
cv <- apply(data4_nonzero,1,sd)/rowMeans(data4_nonzero)# 计算变异系数(标准差除以平均值)度量基因表达变化幅度
data4_use_log2 <- data4_nonzero[order(cv,decreasing = T),]# 根据变异系数排序

# 计算中值绝对偏差 (MAD, median absolute deviation)度量基因表达变化幅度
# 在基因表达中，尽管某些基因很小的变化会导致重要的生物学意义，但是很小的观察值会引入很大的背景噪音，因此也意义不大。
mads <- apply(data4_use_log2, 1, mad)
data4_use_log2 <- data4_use_log2[rev(order(mads)),]
data_var3 <- data4_use_log2[8:10,]#筛选前3列
data_var3_forPCA <- t(data_var3)# 转置矩阵使得每一行为一个样品，每一列为一组变量
dim(data_var3_forPCA)

# 获得样品分组信息
sample <- rownames(data_var3_forPCA)
## One better way to generate group
group <- unlist(lapply(strsplit(sample, "_"), function(x) x[2]))
print(sample[1:7])# 获得样品分组信息
sample <- rownames(data_var3_forPCA)
group <- unlist(lapply(strsplit(sample, "_"), function(x) x[2]))
print(sample[1:7])
print(group[1:7])# 根据分组数目确定颜色变量
colorA <- rainbow(length(unique(group)))# 根据每个样品的分组信息获取对应的颜色变量
colors <- colorA[as.factor(group)]# 根据样品分组信息获得legend的颜色
colorl <- colorA[as.factor(unique(group))]
scatterplot3d(data_var3_forPCA[,1:3], color=colors,pch=16)
legend("topright", legend=levels(as.factor(group)), col=colorl, pch=16, xpd=T, horiz=F, ncol=6)

#2dPCA
data4_use_log2_t <- t(data4_use_log2)
# Add group column for plotting
data4_use_log2_label <- as.data.frame(data4_use_log2_t)
data4_use_log2_label$group <- group
pca <- prcomp(data4_use_log2_t, scale=T)
percentVar <- pca$sdev^2 / sum( pca$sdev^2)
print(str(pca))
autoplot(pca, data=data4_use_log2_label, colour="group") + xlab(paste0("PC1 (", round(percentVar[1]*100), "% variance)")) + ylab(paste0("PC2 (", round(percentVar[2]*100), "% variance)")) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                      panel.grid.minor = element_blank()) + theme(legend.position="right")

#3d pca根据分组数目确定颜色变量
colorA <- rainbow(length(unique(group)))
# 根据每个样品的分组信息获取对应的颜色变量
colors <- colorA[as.factor(group)]
# 根据样品分组信息获得legend的颜色
colorl <- colorA[as.factor(unique(group))]
pc <- as.data.frame(pca$x)
scatterplot3d(x=pc$PC1, y=pc$PC2, z=pc$PC3, pch=16, color=colors, xlab=paste0("PC1 (", round(percentVar[1]*100), "% variance)"),
              ylab=paste0("PC2 (", round(percentVar[2]*100), "% variance)"), zlab=paste0("PC3 (", round(percentVar[3]*100), "% variance)"), angle=45)
legend("topright", legend=levels(as.factor(group)), col=colorl, pch=16, xpd=T, horiz=F, ncol=6)

#hclust
out.dist=dist(t(data4_nonzero),method="euclidean")
out.hclust=hclust(out.dist,method="complete")
plot(out.hclust)

#volcano plot p_p_vs_p_h
setwd("~/yinxiaolan")
cts<-read.table("tumor_vs_normal",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("tumor",7),rep("normal",7)), levels = c("tumor","normal"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("tumor","normal"))
dds$condition <- relevel(dds$condition, ref = "normal")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"tumor_vs_normal_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","tumor","normal"))
resultsNames(dds2)
write.csv(res,"tumor_vs_normal_deg.csv")
res
res <- results(dds2, name="condition_tumor_vs_normal")
res <- results(dds2, contrast=c("condition","tumor","normal"))

#Log fold change shrinkage for visualization and ranking
resultsNames(dds2)
resLFC <- lfcShrink(dds2, coef="condition_tumor_vs_normal", type="apeglm")
resLFC
library("BiocParallel")
register(MulticoreParam(4))

resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds2, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

# (unevaluated code chunk)
library("IHW")
resIHW <- results(dds2, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

DESeq2::plotMA(res,ylim=c(-10,15))
DESeq2::plotMA(resLFC,ylim=c(-10,15))#用DESeq2包中的plotMA
