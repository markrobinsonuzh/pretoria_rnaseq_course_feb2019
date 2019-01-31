# to activate the tools needed:
source activate rnaseq

# In this tutorial we will do:
# 1) exploratory plots;
# 2) DGE analyses via edgeR;
# 3) DGE analyses via DESeq2;
# 4) compare the results via a Venn diagram.

# download the tximport count matrices from github:
https://github.com/markrobinsonuzh/pretoria_rnaseq_course_feb2019/blob/master/tutorials%20Simone/txi_counts_matrices.RData

wget https://github.com/markrobinsonuzh/pretoria_rnaseq_course_feb2019/blob/master/tutorials%20Simone/txi_counts_matrices.RData
wget https://github.com/SimoneTiberi/BG4-2018/blob/master/txi_counts_matrices.RData
wget https://github.com/SimoneTiberi/BG4-2018/blob/master/txi_counts_matrices.RData?raw=true
# raw=true
# clone the Github repo:
get clone https://github.com/SimoneTiberi/BG4-2018.git

##########################################################################################
# 1) Exploratory plots
##########################################################################################
rm(list = ls())
load("txi_counts_matrices.RData")
ls()
str(txi_gene); # gene-level count matrix, abundance and average transcript length

counts = txi_gene$counts
head( counts ) # gene-level count matrix

# 2 groups: 3 females and 3 males
# 1 additional covariate: time (6, 15 and 18).
samples = c("F06", "F15", "F18",
            "M06", "M15", "M18")

# 3-D MDS plot:
library(edgeR); library(rgl); library(limma)

### female vs male:
group_id = c( rep("female", 3), rep("male", 3) )
cols = ifelse(group_id == "female", "red", "blue")

# Creates a DGEList object:
d_gene = DGEList(counts = counts, samples = samples);
d_gene
# calcNormFactors computes normalizing factors to scale w.r.t. the library size:
d_gene = calcNormFactors(d_gene)
d_gene

# log-cpms (not to over-weight highly expressed genes)
CPM = cpm(d_gene, log = TRUE)
colnames(CPM) = samples

# 2D MDS plot based on Gene-level counts:
plotMDS(CPM, labels = samples, col = cols, main = "Gene MDS")
# not a clear clustering structure between F and M.
# instead, a nice separation driven by time: we MUST account for time, when testing btw F and M.

# PCA plot:
library(ggfortify)
autoplot(prcomp(t(CPM)), # CPM
         colour = cols,
         label = TRUE)

# k-means clustering with 2 and 3 groups:
km_2 = kmeans(t(CPM), 2)
km_2$cluster
F06 F15 F18 M06 M15 M18 
1   1   2   1   1   2 
# if we consider a clustering in 2 groups, the separation is driven by time: 18 vs (6 and 15)
# Correlation plots between samples's expression patterns:
library(corrplot)
corrplot(cor(CPM), method="color" )
# Again: not a clear male/femal separation, and 18 have low correlation with the other times.

# heatmap based on the most variable GENES:
var_genes <- apply(CPM, 1, var)
# we select the 10^3 most variable log-cpms
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
# Subset CPM matrix (only keep the most variables log-cpm)
highly_variable_lcpm <- CPM[select_var,]
library( gplots) 
# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", 
          main="Top 1000 variable Genes")
# Again gender does not create a clear separation;
# time=18 seems to be the most distinguishing feature.
# We want to model time as a discrete covariate (15 is "more similar" to 6 than to 18).

# hierarchical clustering
# samples are clustered on the top (columns)
# genes are clustered on the side (rows)
# useful plot to see:
# 1) how samples cluster;
# 2) how genes cluster together, often genes with a similar function or in the same pathway cluster together.

# Note about hierarchical clustering: the length of the tree branches is proportional to the dissimilarity measure.

# scatterplot of pairs of covariates:
pairs(highly_variable_lcpm)
# same story as above.

##########################################################################################
# 2) DGE analysis with edgeR
##########################################################################################
# edgeR vignette: https://www.google.com/search?q=edgeR+vignette&rlz=1C5CHFA_enCH779CH779&oq=edgeR+vignette&aqs=chrome..69i57j69i60l3j35i39j0.2312j0j7&sourceid=chrome&ie=UTF-8

rm(list = ls())
load("txi_counts_matrices.RData")
ls()
str(txi_gene); # gene-level count matrix, abundance and average transcript length

cts     <- txi_gene$counts
normMat <- txi_gene$length
normMat <- normMat/exp(rowMeans(log(normMat)))
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

samples = c("F06", "F15", "F18",
            "M06", "M15", "M18")
group_id = c( rep("female", 3), rep("male", 3) )
y <- DGEList(cts, samples = samples, group = group_id)
y <- scaleOffset(y, t(t(log(normMat)) + o))
y

# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

y <- calcNormFactors(y)
y
# norm.factors column
# a gene and sample specific offset was created to account for the average transcript length.

####################################################################################
# Ignore time, group is the only covariate:
####################################################################################
group_id
design <- model.matrix(~group_id)
design # design matrix of our model, the only covariate is the group (A or B)

# first, we estimate the dispersion
y <- estimateDisp(y, design)
y$common.dispersion
[1] 0.4563105
# individual dispersion estimates:
head(y$tagwise.dispersion)
# moderated dispersion estimates:
head(y$trended.dispersion)

# exact likelihood ratio test test:
# (the exact test is only applicable to experiments with a single factor)
lrt_exact = exactTest(y)
lrt_exact

# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2) # coef specifies the coefficient to test.
lrt

# compare the two results (exact and approximate p.value):
plot(lrt_exact$table$PValue, lrt$table$PValue )
# the exact test is only applicable to experiments with a single factor.

topTags(lrt_exact)
# FDR represents the adjusted p.value (BH)
# to see all results: topTags(lrt, n = Inf)
# We can also sort by fold change (FC)
topTags(lrt_exact, sort.by="logFC")
exp(7.397927)
# 1632.597
# Gene "TRINITY_DN3576_c0_g1" is expressed ~1,600 times more in males than in females.
# it's the same gene as above, time didn't affect this gene.
design

# plot biological coefficient of variation vs log(CPM), where CPM = Counts per million
plotBCV(y)

# plot FC (fold change) vs log(CPM)
plotSmear(y)
# there is more variability for low counts: on the left of the plot, we observe the most extreme FCs.

save(lrt_exact, lrt, file = "Results/Results_edgeR.Rdata")

####################################################################################
# Add time as a covariate:
####################################################################################
samples
time = factor(c(6, 15, 18, 6, 15, 18))
time # we model time as a cathegorical variable with 3 levels (see exploratory plots)
design <- model.matrix(~group_id + time)
design # design matrix of our model, the only covariate is the group (A or B)

# first, we estimate the dispersion
y <- estimateDisp(y, design)
y$common.dispersion
[1] 0.1878733
# individual dispersion estimates:
head(y$tagwise.dispersion)
# moderated dispersion estimates:
head(y$trended.dispersion)

# exact likelihood ratio test test:
# (the exact test is only applicable to experiments with a single factor)
# lrt_exact = exactTest(y)
# lrt_exact
# by default it tests the first 2 coefficients: intercept vs male

# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit <- glmFit(y, design)
design
lrt <- glmLRT(fit, coef=2) # coef specifies the coefficient to test.
lrt

lrt_time18 <- glmLRT(fit, coef=4) # coef specifies the coefficient to test.
lrt_time18 # here we test time18 vs (baseline) time6.

# compare the two results (exact and approximate p.value):
plot(lrt_exact$table$PValue, lrt$table$PValue )
# the exact test is only applicable to experiments with a single factor.

topTags(lrt_exact)
# FDR represents the adjusted p.value (BH)
# to see all results: topTags(lrt, n = Inf)
# We can also sort by fold change (FC)
topTags(lrt_exact, sort.by="logFC")
exp(7.397927)
# 1639.011
# Gene "TRINITY_DN3576_c0_g1" is expressed ~1,600 times more in males than in females.
design

# the p-value, or adjusted p-value (FDR), represents the evidence that there is an effect (statistical significance).
# The p-value is mainly influenced by
# 1) the magnitude of the effect (bigger effects are easier to detect);
# 2) the data available: typically, the more counts, the easier it is to detect an effect.
# However, it does tell us how "big" the effect is: 
# a small p-value indicates a "clear/compelling" differnece, but not that this difference is big.
# The (estimated) FC (or log-FC), instead, measures how strong the change is between conditions (biological relevance).
# FC = 10, indicates that one gene is (estimated to be) 10 times more expressed in one condition than in the other,
# but it may not be significant: we typically consider both FC and p-value (or FDR).

# plot biological coefficient of variation vs log(CPM), where CPM = Counts per million
plotBCV(y)

# plot FC (fold change) vs log(CPM)
plotSmear(y)
# there is more variability for low counts: on the left of the plot, we observe the most extreme FCs.

save(lrt_exact, lrt, file = "Results/Results_edgeR_InclTime.Rdata")

####################################################################################
# Compare results, before and after adding time:
####################################################################################
load("Results/Results_edgeR_InclTime.Rdata")
lrt_exact_time = lrt_exact
lrt_time = lrt

rm(lrt); rm(lrt_exact)
load("Results/Results_edgeR.Rdata")

plot(lrt_time$table$PValue, lrt$table$PValue)
abline(0,1, col = "red", lwd = 3)

dim(lrt_time$table); dim(lrt$table)

sum(lrt_time$table$PValue < 0.05)
sum(lrt$table$PValue < 0.05)

# accounting for time increases the number of significant genes.

####################################################################################
# Differential expression above a fold-change threshold:
####################################################################################
# The classical DGE approach tests if FC = 1 vs |FC| > 1, or analogously, H0: log(FC) = 0 vs H1: |log(FC)| > 1
# But a statistically significant FC of 1.1, might not be biologically relevant.
# We can be more conservative and detect biologically significant FCs, instead, by testing:
# H0: |log(FC)| < theshold vs H1: |log(FC)| > theshold.

# below threshold = lfc
design
fit <- glmQLFit(y, design)
tr <- glmTreat(fit, coef=2, lfc=1)
topTags(tr)

##########################################################################################
# 3) DGE analysis with DESeq2
##########################################################################################
# DESeq2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

rm(list = ls())
load("txi_counts_matrices.RData")
ls()
str(txi_gene); # gene-level count matrix, abundance and average transcript length

library(DESeq2)
# The user should make sure the rownames of sampleTable align with the colnames 
# of txi$counts, if there are colnames. The best practice is to read 
# sampleTable from a CSV file, and to construct files from a column of sampleTable, 
# as was shown in the tximport examples above.

####################################################################################
# Ignore time, group is the only covariate:
####################################################################################
sampleTable <- data.frame(condition = factor(rep(c("female", "male"), each = 3)) )
sampleTable
rownames(sampleTable) <- colnames(txi_gene$counts)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length:
dds <- DESeqDataSetFromTximport(txi_gene, sampleTable, design= ~condition )
# using counts and average transcript lengths from tximport

# filter lowly abundant genes (less than 10 counts)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# BiocParallel::MulticoreParam(workers = 4) parallelizes computations on 4 cores.
# run everything in 1 command:
dds <- DESeq(dds, BPPARAM = BiocParallel::MulticoreParam(workers = 4))
resultsNames(dds) # lists the coefficients

# dispersion estimates, in black and blue, before and after shrinking to the common trend:
plotDispEsts(dds)

# object 'name' needs to match the coefficient above:
res <- results(dds, name="condition_male_vs_female",
               BPPARAM = BiocParallel::MulticoreParam(workers = 4))
res

# MA-plot
plotMA(res, ylim=c(-2,2)) # significant genes are in red.

# Plot counts:
# normalized counts, separated per group, for the most significant gene:
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# normalized counts, separated per group, for the gene with the biggest FC
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="condition")
# a lot of uncertainty in this gene: the p.value is NA (see below)

# order results by adjusted p-value (FDR)
res[ order(res$padj), ]

# order results by adjusted p-value (FDR)
res[ order( abs(res$log2FoldChange), decreasing = TRUE), ]
# log2(FC) = yy, so FC = 2^yy
# 2^7.56815945483854
# 189.7768
# TRINITY_DN3576_c0_g1 is estimated to have ~190 times higher expression in males compared to females.

# or to shrink log fold changes association with condition:
res_shrink <- lfcShrink(dds, coef="condition_male_vs_female", 
                        BPPARAM = BiocParallel::MulticoreParam(workers = 4))
res_shrink[ order( abs(res_shrink$log2FoldChange), decreasing = TRUE), ]
# TRINITY_DN3576_c0_g1 is also top of the list here for FC.
# 2^5.77070461040883
# 54.59529
# FC of TRINITY_DN3576_c0_g1 has changed a lot (in a conservative manner).

save(res, file = "Results/Results_DESeq2.RData")
rm(res)

####################################################################################
# Add time as a covariate:
####################################################################################
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length:
sampleTable <- data.frame(condition = factor(rep(c("female", "male"), each = 3)), 
                          time = factor(rep(c("6", "15", "18"),  2)) )
sampleTable
rownames(sampleTable) <- colnames(txi_gene$counts)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length:
dds <- DESeqDataSetFromTximport(txi_gene, sampleTable, design= ~condition + time)
# using counts and average transcript lengths from tximport


# filter lowly abundant genes (less than 10 counts)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# BiocParallel::MulticoreParam(workers = 4) parallelizes computations on 4 cores.
# run everything in 1 command:
dds <- DESeq(dds, BPPARAM = BiocParallel::MulticoreParam(workers = 4))
resultsNames(dds) # lists the coefficients

# dispersion estimates, in black and blue, before and after shrinking to the common trend:
plotDispEsts(dds)

# object 'name' needs to match the coefficient above:
res <- results(dds, name="condition_male_vs_female",
               BPPARAM = BiocParallel::MulticoreParam(workers = 4))
res

# MA-plot
plotMA(res, ylim=c(-2,2)) # significant genes are in red.

# Plot counts:
# normalized counts, separated per group, for the most significant gene:
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# normalized counts, separated per group, for the gene with the biggest FC
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="condition")
# a lot of uncertainty in this gene: the p.value is NA (see below)

# order results by adjusted p-value (FDR)
res[ order(res$padj), ]
# TRINITY_DN3576_c0_g1 still is the most significant gene.

# order results by adjusted p-value (FDR)
res[ order( abs(res$log2FoldChange), decreasing = TRUE), ]

save(res, file = "Results/Results_DESeq2_InclTime.RData")

##########################################################################################
# 4) Compare edgeR and DESeq2 significant genes (for a speficied threshold)
# with a Venn diagram
##########################################################################################
# First, we need to sort the results by gene name in order to compare the same gene
# Then, we need to set a significance threshold, e.g. 5%, and select the significant and non-significant genes.
rm(list = ls())
load("txi_counts_matrices.RData")
RES = data.frame(gene_id = rownames(txi_gene$counts))

# edgeR:
load("Results/Results_edgeR_InclTime.Rdata")
res_edgeR = topTags(lrt, n = Inf)$table
match_ = match(RES$gene_id, rownames(res_edgeR))
res_edgeR_ordered = res_edgeR[ match_, ]

# DESeq2:
load("Results/Results_DESeq2_InclTime.RData")
res_DESeq2 = res
match_2 = match(RES$gene_id, rownames(res_DESeq2))
res_DESeq2_ordered = res_DESeq2[ match_2, ]

# select a cut-off (e.g., 0.05)
de_05_edgeR  <- res_edgeR_ordered$FDR < .05
de_05_DESeq2 <- res_DESeq2_ordered$padj < .05

# plot the venn diagramm of the significant genes in the two analyses:
res_05 <- cbind(edgeR=de_05_edgeR, DESeq2=de_05_DESeq2)

head(res_05); tail(res_05)

colSums(res_05, na.rm = T)
# total number of significant genes detected by the two models

vennDiagram(res_05)
# Remember this is not a method evaluation!
# The true status of the genes is unknown here:
# we are just comparing the results to see how similarly the methods' final outputs are.

# what do we conclude ?
# There is a nice agreement between the two methods: most of the genes picked are in common between the two
# However, DESeq2 is less conservative: it has many more significant genes.

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