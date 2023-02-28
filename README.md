# WGCNA-with-HERVS

## Read the data and filter
Filter raw gene matrix and raw erv matrix separately (because ERVs have lower counts and too many are lost when gene filter is applied to them).

```
phe <- read.table("~/Documents/TelescopeDiffExp/samples.design.updated.sv1.171subjects.txt",row.names=1,header=T,check.names=FALSE)
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
If you are starting from here, read the file like this:

```
forInput<-read.table("vsd_norm_filtered_sep_forWGCNA.txt",sep="\t",header=TRUE,row.names=1)
vsdinput<-t(forInput)
```
and then pick up frome there

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
```
### pick power for the model
```
#select power by looking at the plot
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

allowWGCNAThreads()
sft <- pickSoftThreshold(vsdinput,
powerVector = power,
networkType = "unsigned",
verbose = 5)
sft.data <- sft$fitIndices
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
geom_point() +
geom_text(nudge_y = 0.1) +
geom_hline(yintercept = 0.90, color = 'red') +
labs(x = 'Power', y = 'Scale free topology model fit, unsigned R^2') +
theme_classic()
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
geom_point() +
geom_text(nudge_y = 0.1) +
labs(x = 'Power', y = 'Mean Connectivity') +
theme_classic()
grid.arrange(a1, a2, nrow = 2)
```
### select the power (in my case 10) and run the function to create the actual modules.  In this example, I am using 'unsigned'
```
soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor
bwnet <- blockwiseModules(vsdinput,
maxBlockSize = 20000,
TOMType = "unsigned",
power = soft_power,
mergeCutHeight = 0.25,
numericLabels = FALSE,
randomSeed = 1234,
verbose = 3)
cor <- temp_cor
module_eigengenes <- bwnet$MEs
#display numbers of genes in each module
table(bwnet$colors)
```
### plot 
```
mergedColors = labels2colors(bwnet$colors)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
c("unmerged", "merged"),
dendroLabels = FALSE,
addGuide = TRUE,
hang= 0.03,
guideHang = 0.05)
```
### save the module colours
```
write.table(as.data.frame(bwnet$colors),"moduleColors.txt",sep="\t")
```
### separate into ervs and genes so that we can count where ervs are
```
modules <- as.data.frame(bwnet$colors)
modules$erv<-NA
modules$erv[!grepl("ENSG",rownames(modules))] <- 1
modules$erv[grepl("ENSG",rownames(modules))] <- 0
head(modules)
colnames(modules)<-c("module","erv")
table(modules$module,modules$erv)
library(tigerstats)
rowPerc(xtabs(~  module + erv, data=modules))

#sanity check, how many ervs and genes there are
sum(modules$erv ==1, na.rm=TRUE)
sum(modules$erv ==0, na.rm=TRUE)

#get colours for ervs of interest
modules[c("6009","4444","5346","943","2716","1412","4174","W-106"),]
#siginificant DE ervs in KCL BB
BBsigErvs<-modules[c("2574","6253","2476","5480","6251","5346","3145","2658","2845","4310",
"K-50","3656","4243","4444","3949","2982","2449","5833","6285","4184","3346",
"2377","3176","4490","849","5658","704","3922","4613","4146"),]
dim(BBsigErvs)
write.table(BBsigErvs,"BBsigERVS_perModule.txt",quote=FALSE,sep="\t")
```

#### This is specific to KCL BB, where I have phenotype info in two tables that need merging
Skip this and make sure you have traits table which is numeric. The last two lines below are the trick to turn the matrix in numeric
```
phe2 <-read.table("BrainBankPheno.csv",sep=",",header=TRUE,row.names = 1)
test<-merge(phe,phe2,by=0,all=TRUE)
dim(test)
head(test)
rownames(test)<-test$Row.names
traits<-test[,c("Status","AgeCat","Age_Onset","disease_duration_days","Age")]

bckTraits<-traits
charData = apply(as.matrix(traits), 2, as.character);
traits = apply(charData, 2, as.numeric)
```
### Module Trait Correlation
```
nSamples <- nrow(vsdinput)
nGenes <- ncol(vsdinput)
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

#if the first column is the actual rownames, run the following to fix it
heatmap.data <- heatmap.data %>%
column_to_rownames(var = 'Row.names')
```
For the plot, look at the heatmap.data to select the columns which are the phenotype data.  in my case, these are columns 17 to 21, and these will
go to the X-axis.  In y axis put all the columns that are module eigengenes
```
CorLevelPlot(heatmap.data,
x = names(heatmap.data)[17:21],
y = names(heatmap.data)[1:16],
col = c("blue1", "skyblue", "white", "pink", "red"))

```

### Get genes in a module
```
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>%
filter(`bwnet$colors` == 'cyan') %>%
rownames()
```

### Module membership meaasure - to get the leading genes that represent a module
```
module.membership.measure <- cor(module_eigengenes, vsdinput, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvals[1:10,1:10]
```

### Genes correlated with phenotype
```
gene.signf.corr <- cor(vsdinput, traits[,'Age'], use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals %>%
as.data.frame() %>%
arrange(V1) %>%
head(25)
```
