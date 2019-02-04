# to activate the tools needed:
source activate rnaseq

# Overview paper about DTU based on transcript estimated counts:
# Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification 
# https://f1000research.com/articles/7-952/v3

# In this tutorial we will do:
# 1) DTU analyses via DRIMSeq;
# 2) DEU analyses via DEXSeq (actually DTU, since we'll use transcript esrtimated counts);
# 3) compare the results via a Venn diagram.
# 4) 2-stage testing via StageR

##########################################################################################
# 1) DRIMSeq
##########################################################################################
# DRIMSeq vignette: http://bioconductor.org/packages/release/bioc/vignettes/DRIMSeq/inst/doc/DRIMSeq.pdf
rm(list = ls())
load("txi_counts_matrices.RData")
ls()

library(DRIMSeq); library(edgeR)
samples = c("F06", "F15", "F18",
            "M06", "M15", "M18")

# match each transcript with the corresponding gene id
head(GeneToTrans)
matching =  match(rownames(txi$counts), GeneToTrans$Tr_id)
head(GeneToTrans[matching, ]);
head(rownames(txi$counts))
gene_id = GeneToTrans$Gene_id[matching ]
colnames(txi$counts) = samples

counts_df = data.frame(txi$counts, gene_id = gene_id, feature_id = rownames(txi$counts))
head(counts_df)

####################################################################################
# Ignore time, group is the only covariate:
####################################################################################
sampleTable <- data.frame(sample_id = samples,
                      group = c( rep("female", 3), rep("male", 3) ))
sampleTable

# Create a dmDSdata object
d <- dmDSdata(counts = counts_df, samples = sampleTable)
d
head(counts(d), 3)
head(samples(d), 3)
dim(counts(d))

# Nr of features (transcripts) per gene
plotData(d)

table(samples(d)$group)
# filter genes with < 10 counts and transcripts with relative abundance < 0.01
d <- dmFilter(d, min_gene_expr = 100, min_feature_expr = 1, min_feature_prop =  0.1)
d
plotData(d)
# Genes with 1 tr only were removed, ok!
# with 25,954 genes and 6 samples
# min_gene_expr	   Minimal gene expression.
# min_feature_expr Minimal feature expression.
# min_feature_prop Minimal proportion for feature expression. This value should be between 0 and 1.

design_full <- model.matrix(~ group, data = samples(d))
design_full

# Dirichlet-Multinomial precision estimate:
# This takes ~5 mins (on 4 cores)
system.time({d <- dmPrecision(d, verbose = 1, genewise_precision = TRUE, 
                 design = design_full,
                 BPPARAM = BiocParallel::MulticoreParam(workers = 4))})
common_precision(d)
head(genewise_precision(d))
plotPrecision(d)

# We fit the model
# ?dmFit
d <- dmFit(d, design = design_full, verbose = 1, BPPARAM = BiocParallel::MulticoreParam(workers = 4))
d

head(proportions(d))
## Get the Dirichlet-Multinomial (DM) regression coefficients (gene-level)
head(coefficients(d))
## Get the Beta-Binomial (BB) regression coefficients (feature-level)
head(coefficients(d), level = "feature")

# We test the genes
# ?dmTest
# name 'coef' as the coefficient in the design matrix:
design(d)
d <- dmTest(d, coef = "groupmale", verbose = 1, BPPARAM = BiocParallel::MulticoreParam(workers = 4))
d

plotPValues(d) # gene-level p.values
plotPValues(d, level = "feature") # tr-level p.values

results_gene	    = results(d, level = "gene")
results_trancript = results(d, level = "feature")

# Save results:
save(results_gene, results_trancript, file = "Results/DRIMSeq_results.RData")

# plot Results for the TOP gene (smallest p-value):
results_gene <- results_gene[order(results_gene$pvalue, decreasing = FALSE), ]
top_gene_id <- results_gene$gene_id[1]
plotProportions(d, gene_id = top_gene_id, group_variable = "group")
plotProportions(d, gene_id = top_gene_id, group_variable = "group",  plot_type = "lineplot" )
# in ribbonplot you can only see the average proportion (not the sample-specific ones):
plotProportions(d, gene_id = top_gene_id, group_variable = "group",  plot_type = "ribbonplot")

head(results_gene)


##########################################################################################
# 2) DEXSeq
##########################################################################################
# We perform DEU on transcript-level counts, with DEXSeq:
# DEXSeq vignette: https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf

rm(list = ls())
load("txi_counts_matrices.RData")
ls()

library(DEXSeq)
samples = c("F06", "F15", "F18",
            "M06", "M15", "M18")

# match each transcript with the corresponding gene id
head(GeneToTrans)
matching =  match(rownames(txi$counts), GeneToTrans$Tr_id)
head(GeneToTrans[matching, ]);
head(rownames(txi$counts))
gene_id = GeneToTrans$Gene_id[matching ]
colnames(txi$counts) = samples

####################################################################################
# Ignore time, group is the only covariate:
####################################################################################
sampleTable = data.frame( row.names = samples,
                          condition = c(rep("F", 3), rep("M", 3)) )
sampleTable

dxd = DEXSeqDataSet(countData = round( txi$counts ),
                    sampleData=sampleTable,
                    design= ~ sample + exon + 
                      condition:exon,
                    featureID = rownames(txi$counts),
                    groupID = gene_id)
dxd
head( counts(dxd), 5 )
head( featureCounts(dxd), 5 )

# normalization:
dxd = estimateSizeFactors( dxd )

# Parallelize future tasks on 4 cores:
library(BiocParallel)
BPPARAM = BiocParallel::MulticoreParam(workers = 4)

# dispersion estimate:
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
plotDispEsts( dxd )

# NB test for differential exon/transcript usage:
dxd = testForDEU( dxd, BPPARAM=BPPARAM)

# exon/transcript fold change between conditions:
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
plotMA( dxd )

res = DEXSeqResults( dxd) # independentFiltering = FALSE )
res
res[ order(res$padj), ] # the most significant genes on top.
res[ order(abs(res$log2fold_M_F), decreasing = T), ] # the highest abs(log2FC) on top

# perGeneQValue computes an adjusted gene-level p.value,
# acconting for the fact that several transcript/exon tests come from the same gene:
qval = perGeneQValue(res)
head(sort(qval))

# Plot of exon/transcript abundance:
plotDEXSeq( res,geneID="TRINITY_DN10922_c1_g1", cex.axis=1.2, cex=1.3, lwd=2,
            expression = FALSE, norCounts = TRUE)

res_gene = data.frame(gene = names(qval), qval)

save(res, res_gene, file="Results/DEXSeq_results.RData")

rm(dxd)

####################################################################################
# Add time as a covariate:
####################################################################################
sampleTable <- data.frame(condition = factor(rep(c("F", "M"), each = 3)), 
                          time = factor(rep(c("6", "15", "18"),  2)) )
rownames(sampleTable) <- colnames(txi$counts)
sampleTable

dxd = DEXSeqDataSet(countData = round( txi$counts ),
                    sampleData=sampleTable,
                    design= ~ sample + exon + 
                      time:exon
                    + condition:exon,
                    featureID = rownames(txi$counts),
                    groupID = gene_id)
dxd
head( counts(dxd), 5 )
head( featureCounts(dxd), 5 )

# normalization:
dxd = estimateSizeFactors( dxd )

# Parallelize future tasks on 4 cores:
library(BiocParallel)
BPPARAM = BiocParallel::MulticoreParam(workers = 4)

# dispersion estimate:
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
plotDispEsts( dxd )

# NB test for differential exon/transcript usage:
dxd = testForDEU( dxd, BPPARAM=BPPARAM)

# exon/transcript fold change between conditions:
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
plotMA( dxd )

res = DEXSeqResults( dxd) # independentFiltering = FALSE )
res
res[ order(res$padj), ] # the most significant genes on top.
res[ order(abs(res$log2fold_M_F), decreasing = T), ] # the highest abs(log2FC) on top

# perGeneQValue computes an adjusted gene-level p.value,
# acconting for the fact that several transcript/exon tests come from the same gene:
qval = perGeneQValue(res)
head(sort(qval))

# Plot of exon/transcript abundance:
plotDEXSeq( res,geneID="TRINITY_DN10922_c1_g1", cex.axis=1.2, cex=1.3, lwd=2,
            expression = FALSE, norCounts = TRUE)

res_gene = data.frame(gene = names(qval), qval)

save(res, res_gene, file="Results/DEXSeq_results_InclTime.RData")

####################################################################################
# Compare results, before and after adding time:
####################################################################################
load("Results/DEXSeq_results_InclTime.RData")
res_gene_time = res_gene

rm(res_gene)
load("Results/DEXSeq_results.RData")

length(res_gene_time$qval); length(res_gene$qval)

plot(res_gene_time$qval, res_gene$qval)
abline(0,1, col = "red", lwd = 3)

dim(lrt_time$table); dim(lrt$table)

sum(res_gene_time$qval < 0.05)
sum(res_gene$qval < 0.05)
# Similarly to the DGE case, we have many more significant genes after accounting for time.


##########################################################################################
# 3) Compare DRIMSeq, DEXSeq as well as edgeR and DESeq2 significant genes (for a speficied threshold)
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

# DRIMSeq:
load("Results/DRIMSeq_results.RData")
res_DRIMSeq =  results_gene
match_3 = match(RES$gene_id, results_gene$gene_id)
res_DRIMSeq_ordered = res_DRIMSeq[ match_3, ]

# DEXSeq:
load("Results/DEXSeq_results.RData")
res_DEXSeq =  res_gene
match_4 = match(RES$gene_id, res_DEXSeq$gene)
res_DEXSeq_ordered = res_DEXSeq[ match_4, ]


# select a cut-off (e.g., 0.05)
de_05_edgeR    <- res_edgeR_ordered$FDR < .05
de_05_DESeq2   <- res_DESeq2_ordered$padj < .05
de_05_DRIMSeq  <- res_DRIMSeq_ordered$adj_pvalue < .05
de_05_DEXSeq   <- res_DEXSeq_ordered$qval < .05

# plot the venn diagramm of the significant genes in the two analyses:
res_05 <- cbind(edgeR=de_05_edgeR, DESeq2=de_05_DESeq2,
                DRIMSeq=de_05_DRIMSeq, DEXSeq=de_05_DEXSeq)

head(res_05); tail(res_05)

colSums(res_05, na.rm = T)
# total number of significant genes detected by the four models
# Remember that edgeR and DESeq2 perform DGE, while DRIMSeq and DEXSeq DS (diff splicing)!

vennDiagram(res_05)
# Remember this is not a method evaluation!
# The true status of the genes is unknown here:
# we are just comparing the results to see how similarly the methods' final outputs are.

# what do we conclude ?
# There is a decent agreement between the DGE methods and the DS tools.
# It's clear that DGE and DS search for different biological phenomena.

##########################################################################################
# 4) 2-stage testing via StageR:
##########################################################################################
# StageR vignette: http://bioconductor.org/packages/release/bioc/vignettes/stageR/inst/doc/stageRVignette.html