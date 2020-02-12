
#Need to set Session>Set Working Directory>To Source File Location
directory <- "."

#Get all the trapnell count files
sampleFiles <- grep("trapnell_counts_",list.files(directory),value=TRUE)

#Extract the condition from the file names
sampleCondition <- sub("trapnell_counts_(C.*)_R.*","\\1",sampleFiles)

#Put all the information in a table
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)


library("DESeq2")
#Load all the data
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#This runs the whole pipeline: normalization, variance estimation, and testing
ddsHTSeq <- DESeq(ddsHTSeq)

#Plots the dispersion re-estimation procedure
plotDispEsts(ddsHTSeq)

#Produces the results into tables to be saved and to make plots
resHTSeq <- results(ddsHTSeq)

#Save result as csv file, sorted by adjusted pvalue
orderedRes<-resHTSeq[order(resHTSeq$padj),]
write.csv(as.data.frame(orderedRes), file="trapnell_C1_VS_C2.DESeq2.csv")

#MA Plot
plotMA(resHTSeq)

#Vulcano Plot. In red, genes with adjusted pvalue of 0.05
plot(resHTSeq$log2FoldChange,-log2(resHTSeq$pvalue))
points(resHTSeq$log2FoldChange[resHTSeq$padj<0.05],-log2(resHTSeq$pvalue[resHTSeq$padj<0.05]),col="red")

#Plot Histogram of p-values
hist(resHTSeq$pvalue, breaks=0:50/50, xlab="p value", main="Histogram of uncorrected p values")


#Does a variant stabilization procedure useful for some plotting and clustering
#NOT to be used for differential analysis
vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)

#Plot PCA
plotPCA(vsd, intgroup=c("condition"))

#Plot Heatmap
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- sampleTable$sampleName
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#Plot normalized values of top 20 differentaially expressed genes...
log2_norm_counts<-log2(counts(ddsHTSeq,normalized=TRUE)+1)
select<-row.names(orderedRes[1:20,])
  
values=log2_norm_counts[select,]

pheatmap(values,
         scale = "none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_row = 8, #by default is 10
         annotation_names_col = FALSE,
         gaps_col = c(3,6),
         display_numbers = TRUE,
         number_format = "%.2f",         
         #filename = "trapnell_C1_VS_C2.topdiff.png",
         height=12,
         width=6
)

