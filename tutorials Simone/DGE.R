
# In this tutorial we will do:
# 1) exploratory plots;
# 2) DGE analysis via edgeR;
# 3) DGE analysis via DESeq2;
# 4) DGE analysis via DESeq;
# 5) compare the results via a Venn diagram.

cd
mkdir ex_2
cd ex_2

R
# in R:
# data in /gdc_home2/groups/BAG18/5_friday/ex_2/data

# How to install a bioconductor package ?
# install.packages("edgeR") # it won't work/

# correct way:
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")

dir = "/gdc_home2/groups/BAG18/5_friday/ex_2/data/"
samples <- c( paste0("condA_",1:3,".count"),paste0("condB_",1:4,".count"))
samples

##########################################################################################
# 1) Exploratory plots
##########################################################################################

library(edgeR)
counts <- readDGE(paste0(dir, samples) )$counts
colnames(counts) = samples
head(counts); tail(counts)

# we remove the bottom lines (not genes) by filtering the names without "FBgn"
counts <- counts[ grep("^FBgn",rownames(counts)), ]
head(counts); tail(counts)
dim(counts)

group <- c(rep("A",3), rep("B", 4))

library(corrplot)
library( gplots) 

pdf("raw_counts.pdf")
# plot counts of samples belonging to the same group
plot(counts[,5], counts[,6], main = "raw counts", xlab = "B2", ylab = "B3")
abline(0,1)
# plot counts of samples belonging to different groups
plot(counts[,5], counts[,1], main = "raw counts", xlab = "B2", ylab = "A1")
abline(0,1)

# correlation between pairs of samples
corrplot(cor(counts), method="color")
dev.off()
# the counts of the samples belonging to the same group are a lot more similar (closer to the main diagonal) 
# w.r.t. the counts of the samples belonging to different groups.
# the correlation plot should show a higher correlation between the samples of the same group compared to that one between samples of diffrent groups.
# in this specific dataset not very evident.

# we create a DGEList element
y <- DGEList(counts=counts, group=group)
y

logcounts <- cpm(y,log=TRUE)
# we use log counts to diminish the importance of high number counts.

dim(logcounts)
# to speed up the process, we select the 500 most variable genes.

var_genes <- apply(logcounts, 1, var)
head(var_genes)

# we select the 500 most variable log-cpm
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix (only keep the 500 most variables log-cpm)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

colnames(highly_variable_lcpm) = c(paste0("A",1:3), paste0("B", 1:4))

pdf("clustering.pdf")
# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")
dev.off()

# samples are clustered on the top (columns)
# genes are clustered on the side (rows)
# useful plot to see 
# 1) how samples cluster: in this case there is not a clear distinction between the two groups, which suggests the two groups are pretty similar.
# A1 and B1 are separated from the other samples (they show the smallest correlation in the correlation plot).
# 2) how genes cluster together, often genes with a similar function or in the same pathway cluster together.
# top genes distinct from the rest: low expression across all samples.

# Note about hierarchical clustering: the length of the tree branches is proportional to the dissimilarity measure.

##########################################################################################
# 2) DGE analysis with edgeR
##########################################################################################

# we create a DGEList element
y <- DGEList(counts=counts, group=group)
y
# we calculate the normalization factors to scale for the library size
y <- calcNormFactors(y)
y
# norm.factors column

design <- model.matrix(~group)
design # design matrix of our model, the only covariate is the group (A or B)

# first, we estimate the dispersion
y <- estimateDisp(y, design)
y$common.dispersion
[1] 0.03801511
# individual dispersion estimates:
$tagwise.dispersion
[1] 0.03041017 0.76230681 0.22885874 0.02265202 0.01817355 ... 15216 more elements ...
# moderated dispersion estimates:
y$trended.dispersion
[1] 0.05504833 0.12301480 0.12428289 0.03553420 0.03029117 ... 15216 more elements ...


# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
# FDR represents the adjusted p.value (BH)
# to see all results: topTags(lrt, n = Inf)
# We can also sort by fold change (FC)
topTags(lrt, sort.by="logFC")

pdf("edgeR.pdf")
# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c(rep("red",3), rep("blue", 4)))
# remember that MDS (and PCA) plots are useful yet partial representations of the data.

# plot biological coefficient of variation vs log(CPM), where CPM = Counts per million
plotBCV(y)

# plot FC (fold change) vs log(CPM)
plotSmear(y)
# there is more variability for low counts: on the left of the plot, we observe the most extreme FCs.
dev.off()

##########################################################################################
# 3) DGE analysis with DESeq2
##########################################################################################

colData<- data.frame(group = c(rep("A", 3), rep("B", 4)), row.names = colnames(counts))

library(DESeq2)
# create a counts matrix object:
dds = DESeqDataSetFromMatrix( countData = counts,
                              colData = colData,
                              design= ~ group)
                              
# entire analysis performed in:
dds <- DESeq(dds)
# estimating size factors (the normalization for library size)
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates (moderation)
# fitting model and testing

# access results:
res_dds <- results(dds)
res_dds
# log2 fold change (MLE): group B vs A 
# Wald test p-value: group B vs A 
# DataFrame with 15221 rows and 6 columns

resultsNames(res_dds)

# Moderation of log2 fold change (FC)
# "In version 1.16 and higher, we have split the moderation of log2 fold changes into a separate function, lfcShrink"
res_dds <- lfcShrink(dds, coef=2)
res_dds

pdf("DESeq2.pdf")
plotDispEsts(dds)

# Plot of normalised mean versus log2 fold change:
plotMA(res_dds)

# Hist of p.values:
hist(res_dds$pvalue, breaks=100, main = "")
dev.off()


##########################################################################################
# 4) DGE analysis with DESeq
##########################################################################################

library(DESeq)
# create a counts matrix object:
DESeq_counts = newCountDataSet( counts, group )

# estimate the library size of each sample:
DESeq_counts = estimateSizeFactors( DESeq_counts )
sizeFactors( DESeq_counts )

# the original counts:
head( counts( DESeq_counts, normalized=FALSE ) )
# the normalized counts:
head( counts( DESeq_counts, normalized=TRUE ) )

# dispersion estimation (and plot):
DESeq_counts = estimateDispersions(DESeq_counts)

# sample specific dispersion estimates:
head( fData(DESeq_counts) )

# Negative-binomial (NB) model likelihood ratio test (LRT) between conditions:
res = nbinomTest( DESeq_counts, "A", "B" )
# the LRT takes a bit
head(res)
# padj, is an adjusted p.value via BH method.

pdf("DESeq.pdf")
plotDispEsts(DESeq_counts)

# Plot of normalised mean versus log2 fold change:
plotMA(res)

# Hist of p.values:
hist(res$pval, breaks=100, main = "")
dev.off()

##########################################################################################
# 5) Compare edgeR, DESeq and DESeq2 significant genes (at a speficied threshold).
##########################################################################################
# First, we need to sort the results by gene name in order to compare the same gene
# Then, we need to set a significance threshold, e.g. 5%, and select the significant and non-significant genes.

# Remember this is not a method evaluation!
# The true status of the genes is unknown here:
# we are just comparing the results to see how similarly the methods' final outputs are.

res_edgeR = topTags(lrt, n = Inf)$table
match_ = match(res$id, rownames(res_edgeR))
res_edgeR_ordered = res_edgeR[ match_, ]

match_2 = match(res$id, rownames(res_dds))
res_dds_ordered = res_dds[ match_2, ]

de_05_edgeR <- res_edgeR_ordered$FDR < .05
de_05_DESeq <- res$padj < .05
de_05_DESeq2 <- res_dds_ordered$padj < .05

# plot the venn diagramm of the significant genes in the two analyses:
res_05 <- cbind(edgeR=de_05_edgeR,DESeq=de_05_DESeq, DESeq2=de_05_DESeq2)

head(res_05)

pdf("VennDiagram.pdf")
vennDiagram(res_05)
dev.off()

# Most of the genes picked are in common between the two methods.

##########################################################################################
# More on the topic: eQTL, expression Quantitative Trait Loci
##########################################################################################
# Similarly to DGE, we can look for eQTL.

# In DGE,  we test if a gene is associated to DE between conditions.
# In eQTL, we test if a gene is associated to DE between phenotypes (defined by SNPs):
# the grouping is defined by the phenotypes, hence gene-level information about SNPs locations is essential.
# We search for the SNPs which define phenotypes associated to changes in the gene expression levels.

# The modelling assumptions (the negative-binomial) and the testing procedure are unchanged.

# In eQTL we apply many more tests than in DGE.
# In DGE we have 1 test per gene;
# in eQTL we have many tests per gene: for every gene we test several SNPs.
# Typically, for every gene, we do not test all possible SNPs.
# Indeed, we constrain the search space (for computational reasons and to diminish the number of tests)
# to the SNPs up to a certain distance to the gene of interest, which are more likely to be associated to the gene of interest.
