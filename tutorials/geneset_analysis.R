
# load packages
library(limma)
library(msigdbr)

# load data
d <- readRDS("edgeR_object.rds")

(design <- d$design)

contrast <- makeContrasts("TRT2-CTRL",
                          levels = colnames(d$design))
contrast


m_df = msigdbr(species = "Mus musculus") %>%
  filter(gs_cat %in% c("H", "C5", "C7"))

sets <- split(m_df$gene_symbol, m_df$gs_name)
n <- sapply(sets, length)
sets <- sets[n >= 20 & n < 1000]
length(sets)


  inds <- ids2indices(sets, d$genes$symbol,
                      remove.empty=TRUE)

  v <- voom(d, design, plot = TRUE)
  f <- lmFit(v, design)
  f <- eBayes(f)
  cf <- contrasts.fit(f, contrast)
  cf <- eBayes(cf)


  cm <- camera(v, inds, design=design,
               contrast=contrast)


tstat <- cf$t[,1]
indices <- inds[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

barcodeplot(tstat, index = indices, cex.main=.8,
              quantiles=c(-1,1)*qt(.95, df=6))

