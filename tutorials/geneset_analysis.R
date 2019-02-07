
## load packages
library(limma)
library(msigdbr)

## load data
d <- readRDS("edgeR_object.rds")

# explore design matrix / contrast / data
(design <- d$design)

(contrast <- makeContrasts("TRT2-CTRL",
                          levels = colnames(d$design)))

d$samples
dim(d)

## do DE analysis with voom
v <- voom(d, design, plot = TRUE)

f <- lmFit(v, design)
f <- eBayes(f) # moderate variance

cf <- contrasts.fit(f, contrast)
cf <- eBayes(cf)

topTable(cf)

# do a sanity check
barplot(v$E["ENSMUSG00000035385.Ccl2",], las=2)

# get gene sets from MSIGDB and take the HALLMARK sets
m_df <- msigdbr(species = "Mus musculus")
keep <- m_df$gs_cat=="H"
m_df <- m_df[keep,]

sets <- split(m_df$gene_symbol, m_df$gs_name)
length(sets)

# match up the gene names from genesets to dataset
inds <- ids2indices(sets, d$genes$symbol,
                    remove.empty=TRUE)

## Run 'camera'
cm <- camera(v, inds, design=design,
             contrast=contrast)
  
  
# pull out t-statistics, make 'barcodeplot'
tstat <- cf$t[,1]
indices <- inds[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

barcodeplot(tstat, index = indices,
            quantiles=c(-1,1)*qt(.95, df=6))

