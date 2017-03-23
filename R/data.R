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
#' typically provided by gene_miRNA_interaction_filter
"mir_interactions"

#' miRNA family mapping table
#'
#' @source ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.gz
#' @format A data frame used to map miR.family ids to MiRBase.Accession ids
"mir_family_info"

#' covariance matrices under the null hypothesis that sensitivity correlation is zero
#'
#' @format A list (different gene-gene correlations k) of lists (different number of miRNAs m) of
#' covariance matrices
"cov.matrices"

#' mircode predicted miRNA gene interactions
#' @source http://www.mircode.org/download.php
#' @format A matrix gene symbols vs miRNA family names. 1 if interaction is predicted, 0 otherwise
"mircode"

#' mircode predicted miRNA gene interactions
#' @source http://www.mircode.org/download.php
#' @description This is a subset of the mircode highly conserved interactions for non-coding RNAs
#' and showing very high total conservation > 20%
#' @format A matrix gene symbols vs miRNA family names. 1 if interaction is predicted, 0 otherwise
"mircode_noncoding"