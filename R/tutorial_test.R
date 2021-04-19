## ---- warning=FALSE, message=FALSE--------------------------------------------
library(SPONGE)

## ---- eval=FALSE--------------------------------------------------------------
#  head(gene_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(gene_expr[1:5,1:8])

## ---- eval=FALSE--------------------------------------------------------------
#  head(mir_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(mir_expr[1:5,1:5])

###TEST BATCH CORRECTION###
#USUALLY WE NEED A LIST OF BATCHES HERE, FOR NOW ITS JUST RANDOM
#gene_expr <- sponge_batch_correction(expression=gene_expr)
#mir_expr <- sponge_batch_correction(expression=mir_expr)

knitr::kable(gene_expr[1:5,1:8])
knitr::kable(mir_expr[1:5,1:5])

###TEST BATCH CORRECTION###

## ---- eval = FALSE------------------------------------------------------------
#  head(targetscan_symbol)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(targetscan_symbol[1:5,1:5])

## ---- warning=FALSE,message=FALSE---------------------------------------------
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
gene_expr = gene_expr,
mir_expr = mir_expr,
mir_predicted_targets = targetscan_symbol)
batches= sample(1:3, nrow(gene_expr), replace=T)
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr,
  mir_expr = mir_expr,
  mir_predicted_targets = targetscan_symbol,
  batches = batches)

## -----------------------------------------------------------------------------
genes_miRNA_candidates[1:2]

## ---- message=FALSE, warning=FALSE--------------------------------------------
ceRNA_interactions <- sponge(gene_expr = gene_expr,
                        mir_expr = mir_expr,
                        mir_interactions = genes_miRNA_candidates)

## ---- eval=FALSE--------------------------------------------------------------
#  head(ceRNA_interactions)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(ceRNA_interactions))

## ---- message=FALSE, warning=FALSE--------------------------------------------
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, number_of_samples = nrow(gene_expr))

## ---- fig.width = 12, fig.height = 7------------------------------------------
sponge_plot_simulation_results(mscor_null_model)

## ---- message=FALSE, error=FALSE----------------------------------------------
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions,
                                                   null_model = mscor_null_model)

## ---- eval=FALSE--------------------------------------------------------------
#  head(ceRNA_interactions_sign)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(ceRNA_interactions_sign))

## -----------------------------------------------------------------------------
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < 0.2),]

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(ceRNA_interactions_fdr))

## -----------------------------------------------------------------------------
sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)

## -----------------------------------------------------------------------------
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)


## -----------------------------------------------------------------------------
ceRNA_interactions_fdr_weight <- ceRNA_interactions_fdr
ceRNA_interactions_fdr_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)

## ---- fig.height = 7, fig.width = 7, warning=FALSE,message=FALSE,error=FALSE----
sponge_plot_network_centralities(weighted_network_centralities, top = 1)

## ---- fig.height = 7, fig.width = 7, warning=FALSE,message=FALSE,error=FALSE----
sponge_plot_network_centralities(weighted_network_centralities, measure = "btw", top = 1)

## ---- eval=FALSE--------------------------------------------------------------
#  library(doParallel)
#  library(foreach)
#
#  num.of.cores <- 2 #many more on a compute cluster
#
#  #if you want to use logging
#  #logging.file <- "where_my_log_file_should_go.log"
#  logging.file <- NULL
#
#  cl <- makeCluster(num.of.cores, outfile=logging.file)
#  registerDoParallel(cl)
#
#  genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
#  gene_expr = gene_expr,
#  mir_expr = mir_expr,
#  mir_predicted_targets = targetscan_symbol)
#
#  ceRNA_interactions <- sponge(
#  gene_expr = gene_expr,
#  mir_expr = mir_expr,
#  mir_interactions = genes_miRNA_candidates)
#
#  stopCluster(cl)
#

## ---- eval = FALSE------------------------------------------------------------
#  library("BiocParallel")
#  register(DoparParam(), default = TRUE)

## ---- message=FALSE, warning=FALSE, error=FALSE-------------------------------
more_covariance_matrices <- sample_zero_mscor_cov(m = 1,
                      number_of_solutions = 10,
                      gene_gene_correlation = 0.5)

## ---- message=FALSE, warning=FALSE, error=FALSE, fig.width = 7, fig.height = 7----
mscor_coefficients <- sample_zero_mscor_data(cov_matrices = more_covariance_matrices,
number_of_samples = 200, number_of_datasets = 100)

hist(unlist(mscor_coefficients), main = "Random distribution of mscor coefficients")