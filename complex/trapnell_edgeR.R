library(edgeR)
library(ggplot2)
library(pheatmap)

countdata<-read.table("trapnell_counts.tab",sep="\t",header=T)
row.names(countdata) <- countdata$Gene
countdata$Gene<-NULL
mygroups <- c("C1","C1","C1","C2","C2","C2")

#Simple, no filtering
y <- DGEList(counts=countdata[,1:6], genes=rownames(countdata), group = mygroups)
y <- calcNormFactors(y)
y$samples
mds <- plotMDS(y, col=c(rep("black",3), rep("red",3)))

#Alternate pca
pc<-prcomp(countdata, retx=TRUE)
percentages<-(pc$sdev)^2 / sum(pc$sdev^2) 
pc<-pc$rotation
pc <- cbind.data.frame(pc,mygroups)
qplot(PC1, PC2, data = as.data.frame(pc), color=mygroups)+geom_point(size=4)+theme(plot.title = element_text(hjust = 0.5, size=20))+labs(x=paste('PC1 (',round(percentages[1], digits=4)*100,'%)'), y=paste('PC2 (',round(percentages[2], digits=4)*100,'%)'), title="PCA")
dev.copy(pdf, "PCA.pdf", width=9, height=7)
dev.off()

#Top genes
y <- estimateDisp(y)
et <- exactTest(y)
topgenes<-topTags(et,n=dim(countdata)[[1]])
plot(topgenes$table$logFC, -log10(topgenes$table$FDR), col=ifelse(topgenes$table$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
hist(topgenes$table$PValue,breaks=20,xlab="P Value",col="royalblue",ylab="Frequency",main="P-value distribution")

normcnt <- round(cpm(y, normalized.lib.sizes=T))
rownames(normcnt)<-rownames(countdata)
topgenes$table <- merge(normcnt, topgenes$table, by.x="row.names", by.y="genes")
names(topgenes$table)[1] <- "Gene.ID"
topgenes$table <- topgenes$table[order(topgenes$table$FDR),]

write.table(topgenes$table, "edgeR_analysis.trapnell_normalized_counts.tab", sep="\t", quote = F, row.names = F)

topgenes_50<-topTags(et, n=50)
#different types of sorting for output
#topgenes_50$table <- topgenes_50$table[order(topgenes_50$table$logFC),]
#topgenes_50$table <- topgenes_50$table[order(topgenes_50$table$logCPM),]
topgenes_50$table <- topgenes_50$table[order(topgenes_50$table$FDR),]
pheatmap(normcnt[as.integer(rownames(topgenes_50$table)),], scale="none", cluster_rows = F, cluster_cols = F, legend = T, main="Top 50 DE Genes")+theme(plot.title = element_text(hjust = 0.9, size=10))
dev.copy(pdf, "heatmap.pdf", width=3, height=10)
dev.off()


for (i in 1:3){
  print(topgenes[i,])
  gene<-as.data.frame(topgenes[i,][,2:7])
  gene_name<-as.data.frame(topgenes[i,][,1])
  print(gene_name)

  stripchart(gene, vertical=TRUE, main=gene_name, pch=16, cex=1.4, col='royalblue', ylab = "Counts per Million")
  dev.copy(pdf, paste(gene_name,".pdf"), width=9, height=7)
  dev.off()
}
  
