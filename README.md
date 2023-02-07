# WGCNA-with-HERVS

## Read the data and filter
Filter raw gene matrix and raw erv matrix separately (because ERVs have lower counts and too many are lost when gene filter is applied to them).

```
gene_cts <- read.table("~/Documents/TelescopeDiffExp/geneMatrixEdited.txt",row.names=1,header=T,check.names=FALSE)
gene_cts<-gene_cts[,rownames(phe)]
cts <- read.table("~/Documents/BrainBank_DE_HERVS/merged_ervs_raw.txt",row.names=1,header=T,check.names=FALSE)
cts<-cts[,rownames(phe)]

phe$Sex<- as.factor(phe$Sex)
phe$Status<- as.factor(phe$Status)
phe$AgeCat<- as.factor(phe$AgeCat)

dds_genes <- DESeqDataSetFromMatrix(countData = gene_cts, colData = phe , design = ~ Sex + AgeCat + Status)
dds_genes

#filter: count of more than 20 in more than 75% of people
idx <- rowSums( counts(dds_genes, normalized=F) >= 20 ) >= 128
dds_genes <- dds_genes[idx,]
dds_genes

dds_ervs <- DESeqDataSetFromMatrix(countData = cts, colData = phe , design = ~ Sex + AgeCat + Status)
dds_ervs

#filter: count of more than 5 in more than 20 of people
idx <- rowSums( counts(dds_ervs, normalized=F) >= 5 ) >= 20
dds_ervs <- dds_ervs[idx,]
dds_ervs

#merge thesee two filtered matrices

#make sure that the column names are identical
identical(colnames(counts(dds_genes)),colnames(counts(dds_ervs)))
#[1] TRUE
gene_erv_cts<-rbind(counts(dds_genes),counts(dds_ervs))
rm(dds_genes)
rm(dds_ervs)
```

## vsd normalise and save the normalised table for WGCNA
```
dds_genes <- DESeqDataSetFromMatrix(countData = gene_erv_cts, colData = phe , design = ~ Sex + AgeCat + Status)
vsd <- vst(dds_genes, blind = TRUE)
dim(assay(vsd))
#[1] 19044   171
write.table(assay(vsd),file="vsd_norm_filtered_sep_forWGCNA.txt",sep="\t",quote=FALSE)
```

## WGCNA
Load the libraries
```
library(tidyverse)
library(magrittr)
library(WGCNA)
library(CorLevelPlot)
```
Either continue with the assay(vsd) from above, or read it in from the file we saved earlier

```
vsdinput = t(assay(vsd))

gsg<-goodSamplesGenes(vsdinput)
summary(gsg)
#sanity check, all samples should pass this check

plot(hclust(dist(vsdinput),method="average"))
