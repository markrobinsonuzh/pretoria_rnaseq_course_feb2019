# In this tutorial we will compare the performance of three methods in a simulation study, where the "true" status of genes (differential or not) is known.

cd
mkdir ex_4
cd ex_4

R
# in R:

##########################################################################################
# We perform a comparison of 3 methods for differential analyses.
##########################################################################################

load("/gdc_home2/groups/BAG18/5_friday/ex_4/results_for_ex_4.RData")
ls()
# res is a table with the adjusted p.values of 3 methods and the true status of 9840 genes.
head(res); dim(res)
# truth = 0: non differential gene;
# truth = 1: differential gene.

# When assessing the performance of a method, we need to know the "truth"!
# We need to know what genes are differentially expressed or show differential splicing,
# to compare the result of the methods with the "truth".
# This is why methods are tested in simulations, where we know the true status of genes.

# we separate the adjusted p.values and the truth in 2 data frames:
padj_icobra = data.frame( res[,1:3] )
rownames(padj_icobra) = 1:nrow(res)

truth_cobra = data.frame(status = res[,4] )
rownames(truth_cobra) = 1:nrow(res)

library(iCOBRA)
# we create an icobra object
cobra = COBRAData(padj = padj_icobra,
                  truth = truth_cobra,
                  object_to_extend = NULL)

####################################################################################################################################
# icobra plots:
cobraperf <- calculate_performance(cobra, binary_truth = "status")
slotNames(cobraperf)

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
facetted = TRUE)

pdf("iCOBRA.pdf")
# Plot observed true positive rate (TPR) for given adjusted p-value thresholds:
# by default 0.01, 0.05 and 0.1
plot_tpr(cobraplot)

# Plot receiver operating characteristics (ROC) curves.
plot_roc(cobraplot)

# MOST USEFUL PLOTS: plot_fpc and plot_fdrtprcurve
# Plot false positive curves, indicating the number of false positives
# among the top-ranked N variables, for varying values of N.
plot_fpc(cobraplot, stripsize = 1, maxnfdc = 1000, linewidth = 0.5)

# Plot observed true positive rate (TPR) vs observed false discovery rate (FDR), 
# for given adjusted p-value thresholds and/or as curves traced out by considering all threshold values.
plot_fdrtprcurve(cobraplot, linewidth = 1)
# the circles indicate the significant genes at the 1, 5 and 10 % thresholds.

# Plot a Venn diagram showing the overlaps among sets of significant feature 
# for a given adjusted p-value threshold
plot_overlap(cobraplot)

# Precision, recall plot:
library(ROCR)
preds = apply(cobra@padj, 2, function(u) 
  prediction(-u, res[,4]))
perfs <- lapply(preds, performance, "tpr", "fpr")
precisions <- lapply(preds, performance, "prec", "rec")

# Recall vs Precision plot.
par( mfrow = c(1,1) )
plot(precisions[[1]], ylim = c(0,1), xlim = c(0,1))
for(i in 2:length(precisions)){
  plot(precisions[[i]], add = T, col = i)
}
legend(legend = colnames(cobra@padj), 
       lwd = 2, col = 1:ncol(cobra@padj), x = 0, y = 0.6)

dev.off()

# The Precision is defined as the number of true positives 
# over the number of all positives (true positives + false positives).
# Precision = TP/(TP + FP)
# The Recall is defined as the number of true positives 
# over the number of significant genes (true positives + false negatives).
# Recall = TP/(TP + FN)
