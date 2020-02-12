rawdata <- read.delim("edgeR_example2_Cumbie.tab", check.names=FALSE, stringsAsFactors=FALSE,row.names=1)
library(edgeR)
Treat <- factor(substring(colnames(rawdata),1,4))
Treat <- relevel(Treat, ref="mock")
Time <- factor(substring(colnames(rawdata),5,5))

y <- DGEList(counts=rawdata, group=Treat)
keep <- rowSums(cpm(y)>2) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
plotMDS(y)
design <- model.matrix(~Time+Treat)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)
topTags(qlf)
top <- rownames(topTags(qlf))
cpm(y)[top,]
summary(dt <- decideTestsDGE(qlf))
