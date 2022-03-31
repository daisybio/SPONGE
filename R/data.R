#' Gene expression test data set
#'
#' @format A data frame of expression values
#' with samples in columns and genes in rows
"gene_expr"

#' miRNA expression test data set
#'
#' @format A data frame of expression values
#' with samples in columns and miRNA in rows
"mir_expr"

#' miRNA / gene interactions
#'
#' @format A data frame of regression coefficients
#' typically provided by sponge_gene_miRNA_interaction_filter
"mir_interactions"

#' ceRNA interactions
#'
#' @format A data table of ceRNA interactions
#' typically provided by sponge
"ceRNA_interactions"

#' covariance matrices under the null hypothesis that sensitivity correlation
#' is zero
#'
#' @format A list (different gene-gene correlations k) of lists
#' (different number of miRNAs m) of
#' covariance matrices
"precomputed_cov_matrices"

#' A null model for testing purposes
#'
#' @format A list (different gene-gene correlations k) of lists
#' (different number of miRNAs m) of
#' sampled mscor values (100 each, computed from 100 samples)
"precomputed_null_model"

#' mircode predicted miRNA gene interactions
#' @source http://www.mircode.org/download.php
#' @format A matrix gene symbols vs miRNA family names. >=1 if interaction is
#' predicted, 0 otherwise
"mircode_symbol"

#' mircode predicted miRNA gene interactions
#' @source http://www.mircode.org/download.php
#' @format A matrix gene ensembl ids vs miRNA family names. >=1 if interaction
#' is predicted, 0 otherwise
"mircode_ensg"

#' targetscan predicted miRNA gene interactions
#' @source http://www.targetscan.org/vert_71/
#' @format A matrix gene symbols vs miRNA family names. >=1 if interaction
#' is predicted, 0 otherwise
"targetscan_symbol"

#' targetscan predicted miRNA gene interactions
#' @source http://www.targetscan.org/vert_71/
#' @format A matrix gene ensembl ids vs miRNA family names. >=1 if interaction
#' is predicted, 0 otherwise
"targetscan_ensg"

#' example training expression data for spongEffects
#' @format a matrix with gene expression data
"train_cancer_gene_expr"

#' example training miRNA data for spongEffects
#' @format a matrix with miRNA expression data
"train_cancer_mir_expr"

#' example training sample meta data for spongEffects
#' @format a data frame with sample meta data, SUBTYPE must be inside your dataframe
"train_cancer_metadata"

#' example test expression data for spongEffects
#' @format a matrix with gene expression data
"test_cancer_gene_expr"

#' example test miRNA data for spongEffects
#' @format a matrix with miRNA expression data
"test_cancer_mir_expr"

#' example test sample meta data for spongEffects
#' @format a data frame with sample meta data, SUBTYPE must be inside your dataframe
"test_cancer_metadata"

#' example train ceRNA interactions for spongEffects
#' @format (obtained by SPONGE method)
"train_ceRNA_interactions"

#' example train network centralities for spongEffects
#' @format (obtained by SPONGE method)
"train_network_centralities"

#' example potential central nodes
#' @format (downloaded via biomaRt)
"Ensembl_IDs"
