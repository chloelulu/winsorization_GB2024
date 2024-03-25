args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
} else{
  i = as.integer(args[2])
  k = as.integer(args[1])
}


winsor.fun <- function(Y, X, quan) {
  N <- estimateSizeFactors(DESeqDataSetFromMatrix(Y, DataFrame(X),~X))$sizeFactor
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

## We used the code from the original manuscript: https://github.com/xihuimeijing/DEGs_Analysis_FDR/blob/main/scripts/DEGs.R
DESeq2.func <- function(readCount, conditions, fdrcutoff=0.05, part, label = '',outdir){
  suppressMessages(library("DESeq2"))
  ddsCount <- DESeqDataSetFromMatrix(readCount, DataFrame(conditions),~conditions)
  dds<-DESeq(ddsCount)
  res<-results(dds)
  res<-na.omit(res)
  res<-res[res$padj<fdrcutoff,]
  out <- paste0(outdir, part,".DESeq",label,".rst.tsv")
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
  return(res)
}

edgeR.func <- function(readCount, conditions, fdrcutoff=0.05,  part, label = '',outdir){
  suppressMessages(library("edgeR"))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~conditions)
  y <- estimateDisp(y,design) # Considering robust=T when you think your data has potential outlier issue.
  #perform quasi-likelihood F-tests:
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  res<-topTags(qlf, n=nrow(readCount), adjust.method = "BH", sort.by = "PValue", p.value = fdrcutoff)
  out <- paste0(outdir,part,".edgeR",label,".rst.tsv")
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
  return(res)
}

Wilcox.func <- function(readCount, conditions, fdrcutoff=0.05, part, label = '', outdir){
  suppressMessages(library('matrixTests'))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method="TMM")
  count_norm<-cpm(y)
  count_norm<-as.data.frame(count_norm)
  conditionsLevel<-levels(conditions)
  dataMem1<-count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataMem2<-count_norm[,c(which(conditions==conditionsLevel[2]))]
  pvalue<-row_wilcoxon_twosample(dataMem1,dataMem2)$pvalue
  fdr<-p.adjust(pvalue,method = "BH")
  outputRst<-na.omit(as.data.frame(cbind(row.names(count_norm)[which(fdr<fdrcutoff)],pvalue[which(fdr<fdrcutoff)],fdr[which(fdr<fdrcutoff)]),stringsAsFactors=F))
  out <- paste0(outdir,part,".Wilcox",label,".rst.tsv")
  write.table(outputRst, file=out,sep="\t", quote=F,row.names = F,col.names = c("Gene","p-value","FDR"))
  return(res = outputRst)
}

# Datasets were downloaded from original paper's zendo: https://zenodo.org/records/6326786
setwd('/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/OriginalCodeData/permuted_datasets/')
expFiles <- list.files(pattern = 'Count.tsv$') # count table
conFiles.permute <- list.files(pattern = 'ConditionLabel.1000.tsv$') # permuted condition labels

## permuted datasets: no-winsorization
conFile.permute <- conFiles.permute[k]
conditions<-read.table(file=conFile.permute, header = F)
conditions<-factor(t(conditions[i,,drop=F]))

outdir1 = "/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result1/"

DESeq2.func(readCount, conditions, part = paste0(filename,'-permute.',i),outdir = outdir1)
edgeR.func(readCount, conditions, part = paste0(filename,'-permute.',i),outdir = outdir1)
Wilcox.func(readCount, conditions, part = paste0(filename,'-permute.',i), outdir = outdir1)


## permuted datasets: winsorization
readCount.win <- winsor.fun(Y = readCount, X = conditions, quan = 0.93)
DESeq2.func(readCount = readCount.win, conditions, part = paste0(filename,'-permute.qt93.RLE.',i), label = '_winsor',outdir = outdir1)
edgeR.func(readCount = readCount.win, conditions, part = paste0(filename,'-permute.qt93.RLE.',i), label = '_winsor',outdir = outdir1)

readCount.win <- winsor.fun(Y = readCount, X = conditions, quan = 0.95)
DESeq2.func(readCount = readCount.win, conditions, part = paste0(filename,'-permute.qt95.RLE.',i), label = '_winsor',outdir = outdir1)
edgeR.func(readCount = readCount.win, conditions, part = paste0(filename,'-permute.qt95.RLE.',i), label = '_winsor',outdir = outdir1)

readCount.win <- winsor.fun(Y = readCount, X = conditions, quan = 0.97)
DESeq2.func(readCount = readCount.win, conditions, part = paste0(filename,'-permute.qt97.RLE.',i), label = '_winsor',outdir = outdir1)
edgeR.func(readCount = readCount.win, conditions, part = paste0(filename,'-permute.qt97.RLE.',i), label = '_winsor',outdir = outdir1)








