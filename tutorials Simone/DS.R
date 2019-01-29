# In this tutorial we will do:
# 1) DTU via DRIMSeq;
# 2) DEU via DEXSeq;

cd
mkdir ex_3
cd ex_3

R
# in R:

##########################################################################################
# 1) DTU via DRIMSeq
##########################################################################################

library(DRIMSeq)
# package for differential transcript usage (DTU) and splicing quantitative trait loci (sQTL)
# it inputs estimated transcript level counts, obtained via salmon, kallisto etc...
# it assumes a Dirichlet-Multinomial hierarchical structure for the transcript level counts across samples of each group.
# For each gene, it performs a likelihood ratio test (LRT) testing if the average transcript proportions vary between two groups or not.
# The dispersion parameters is kept fixed between the two conditions and is not part of the test: interested lies in the proportions.

# DRIMSeq tests for DTU both at the gene and transcript level.

# DRIMSeq also allows the inclusion of covariates if using a Wald test instead of the LRT.
# Covariates can be of interest, i.e. you might want to test if a covariate is involved in DTU,
# or simple confounding factors, i.e. you might want to remove the effect of a covariate when testing between 2 groups, 
# especially if the covariate is not homogeneous between the 2 groups, e.g. age, if one group is significantly younger than the other.

# Example below taken from DRIMSeq vignette:
library(PasillaTranscriptExpr)
data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")

## Load metadata
pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), 
  header = TRUE, as.is = TRUE)

## Load counts
pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), 
  header = TRUE, as.is = TRUE)
head(pasilla_counts)
# feature_id = trancript_id in this case (it can be SNP_id in sQTL analyses)

## order the data in a data.frame
pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
  								group = pasilla_metadata$condition)
levels(pasilla_samples$group)

# create a dmDSdata object
d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
d 
# d gives a summary of the data: 14112 genes
head(counts(d), 3)
head(samples(d), 3)

## we subset the data to speed up calculations.
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
d <- d[names(d) %in% gene_id_subset, ]
d
# 42 genes

pdf("DRIMSeq.pdf")
plotData(d)

## We filter genes with low count numbers (check parameters and defaults)
# independent filtering step.
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
  min_gene_expr = 10, min_feature_expr = 10)
d
# 26 genes.

## We create the design matrix corresponding to the groups
design_full <- model.matrix(~ group, data = samples(d))
design_full

## We set a seed to make the analysis reproducible
set.seed(123)
## Calculate precision
d <- dmPrecision(d, design = design_full)
d
# it takes a long time if computed on all genes.
common_precision(d)
head(genewise_precision(d))
# it computes a common precision for all genes and a gene-specific precision.

plotPrecision(d)

## We fit the model to the data and maximise the likelihood.
d <- dmFit(d, design = design_full, verbose = 1)
d

## Get fitted proportions
head(proportions(d))
## Get the Dirichlet-Multinomial (DM) regression coefficients (gene-level)
head(coefficients(d))
## Get the Beta-Binomial (BB) regression coefficients (feature-level)
head(coefficients(d), level = "feature")

## We perform a LRT between the two groups, coef indicates the coefficient to test
d <- dmTest(d, coef = "groupKD", verbose = 1)
design(d)
head(results(d), 3)
# adj_pvalue is the adjusted p.value to control the FDR with the BH correction method.

## We access the results via the results function
results_gene	  = results(d, level = "gene")
results_trancript = results(d, level = "feature")
dim(results_gene); dim(results_trancript)
# [1] 26  5
# [1] 114   6
head(results_gene)
head(results_trancript)

# many more transcript level tests: for every gene there are > 1 transcripts.

## Plot results:
plotPValues(d) # gene level
plotPValues(d, level = "feature") # transcript level

## Plot top DTU gene:
# we sort according to the p.value (or adjusted p.value: remember that when adjusting p.values, the order remains unchanged!)
results_gene <- results_gene[order(results_gene$pvalue, decreasing = FALSE), ]
top_gene_id <- results_gene$gene_id[1]

plotProportions(d, gene_id = top_gene_id, group_variable = "group")
plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
  plot_type = "lineplot")
plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
  plot_type = "ribbonplot")

dev.off()

##########################################################################################
# 2) DEU via DEXSeq
##########################################################################################
# DEXSeq performs differential exon usage based on the exon level counts.
# Drawback:  it takes DEU as a surrogate for DTU (they are not the same).
# Advantage: there is no uncertainty about the counts at the exon level.
# It tests for DEU every exon bin.
# A p.value at the gene level is obtained as the minimum p.value of the exon tests (statistically speaking, not very elegant...).
# It is a very popular package and works well.

# Example below taken from DEXSeq vignette:

library(pasilla); library(DEXSeq)

## location of data to load (in the package pasilla):
inDir = system.file("extdata", package="pasilla")
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
basename(countFiles)
## location of the gff file in the pasilla package:
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

## description of the data
sampleTable = data.frame(
   row.names = c( "treated1", "treated2", "treated3", 
      "untreated1", "untreated2", "untreated3", "untreated4" ),
   condition = c("knockdown", "knockdown", "knockdown",  
      "control", "control", "control", "control" ),
   libType = c( "single-end", "paired-end", "paired-end", 
      "single-end", "single-end", "paired-end", "paired-end" ) )

sampleTable

## load data into DEXSeq:
dxd = DEXSeqDataSetFromHTSeq(
   countFiles,
   sampleData=sampleTable,
   design= ~ sample + exon + condition:exon,
   flattenedfile=flattenedFile )
# it takes a bit.

## like for DRIMSeq, we subset the data (for computational reasons) and focus on a subset of genes only:
genesForSubset = read.table( 
  file.path(inDir, "geneIDsinsubset.txt"), 
  stringsAsFactors=FALSE)[[1]]
length(genesForSubset)
# 46 genes selected.

dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]
dxd

head( counts(dxd), 5 )
# 1 row per exon bin: the counts are at the exon bin level instead of the transcript level (as for DTU).

## Estimate size factors (normalization)
dxd = estimateSizeFactors( dxd )
# "Different samples might be sequenced with different depths. In order to adjust for such
# coverage biases, we estimate size factors, which measure relative sequencing depth. DEXSeq
# uses the same method as DESeq and DESeq2, which is provided in the function estimate SizeFactors."

## Estimate dispersions (moderated towards a common trend).
dxd = estimateDispersions( dxd )
# gene-level dispersions are shrinked towards a common trend.

pdf("DEXSeq.pdf")
plotDispEsts( dxd )

## We test for DEU at the exon level.
dxd = testForDEU( dxd )
# gene level p.value is obtained as the minimum p.value among the exons of the gene...not a nice aggregation.

## We estimate the fold change of every exon.
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

## Examine the results:
dxr1 = DEXSeqResults( dxd )
dxr1
# 1 row per exon.

## ----tallyExons------------------------------------------------------------
table ( dxr1$padj < 0.05 )
# nr of significant and non-significant exons.

## Plot Mean expression versus log_2 fold change plot.
## Significant hits (at \\Robject{padj}<0.1) are coloured in red.
plotMA( dxr1, cex=0.8 )
dev.off()


##########################################################################################
# More on the topic: sQTL, splicing Quantitative Trait Loci
##########################################################################################
# Similarly to DTU/DTE/DEU, we can look for sQTL.

# In DTU/DTE/DEU,  we test if a gene is associated to differential splicing between conditions.
# In sQTL, 	       we test if a gene is associated to differential splicing between phenotypes (defined by SNPs):
# the grouping is defined by the phenotypes, hence gene-level information about SNPs locations is essential.
# We search for the SNPs which define phenotypes associated to differential splicing.

# The modelling assumptions (the negative-binomial) and the testing procedure are unchanged.

# In sQTL we apply many more tests than in DTU/DTE/DEU.
# In DTU/DTE/DEU we have 1 test per gene;
# in sQTL we have many tests per gene: for every gene we test several SNPs.
# Typically, for every gene, we do not test all possible SNPs.
# Indeed, we constrain the search space (for computational reasons and to diminish the number of tests)
# to the SNPs up to a certain distance to the gene of interest, which are more likely to be associated to the gene of interest.
