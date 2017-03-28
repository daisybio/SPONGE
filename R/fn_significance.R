compute_simulated_data <- function(cov.matrices,
                                   number.of.datasets = 1e5,
                                   number.of.samples){

    #to reach the necessary number of datasets we need to find out how many
    #datasets to construct from each covariance matrix we have
    number.of.datasets.per.matrix <- ceiling(number.of.datasets / length(cov.matrices))

    scor <- sample(
        unlist(sample_zero_scor_data(cov.matrices = cov.matrices,
                                         number.of.datasets = number.of.datasets.per.matrix,
                                         number.of.samples = number.of.samples)),
        number.of.datasets)

    test_data_dt <- data.table(scor)
    setkey(test_data_dt, scor)

    return(test_data_dt)
}

compute_p_values <- function(partition,
                             cov.matrices = NULL,
                             simulated.data = NULL,
                             number.of.datasets = 1e5,
                             number.of.samples = 100){

    if(nrow(partition) == 1 && is.na(partition$geneA))  return(NULL)
    if(is.null(cov.matrices) & is.null(simulated.data))
        stop("either covariance matrices or simulated data have to be provided.")

    if(!("scor" %in% colnames(partition))) stop("sensitivity correlation missing")

    #check which k and m
    k <- as.character(partition[1,cor_cut])
    m <- as.character(partition[1,df_cut])

    #simulate data using the appropriate covariance matrices
    if(is.null(simulated.data)){
        cov.matrices.partition <- cov.matrices[[m]][[k]]

        test_data_dt <- compute_simulated_data(cov.matrices = cov.matrices.partition,
                                               number.of.datasets = number.of.datasets,
                                               number.of.samples = number.of.samples)
    }
    else{
        test_data_dt <- simulated.data[[m]][[k]]
    }

    number.of.datasets.on.right.side <- length(test_data_dt$scor)

    partition$p.val <- (number.of.datasets.on.right.side -
        test_data_dt[J(partition$scor),
                     .I,
                     roll = "nearest",
                     by = .EACHI]$I) / number.of.datasets.on.right.side

    return(partition)
}

ks <- seq(0.2, 0.90, 0.1)
ms <- seq(1, 8, 1)

#' Compute p-values for SPONGE interactions
#'
#' @param sponge_result A data frame from a sponge call
#' @param number.of.samples number of samples in the expression matrices
#' @param number.of.datasets number of datasets to simulate.
#' @param simulated.data optional, pre-computed simulated data
#' affects the mimimum p-value that can be achieved
#' @import data.table
#' @import gRbase
#' @import MASS
#' @import ppcor
#' @import foreach
#' @import logging
#' @import iterators
#' @seealso sponge_build_null_model
#'
#' @return A data frame with sponge results, now including p-values
#' and adjusted p-value
#' @description This method uses pre-computed covariance matrices that were
#' created for various gene-gene correlations (0.2 to 0.9 in steps of 0.1)
#' and number of miRNAs (between 1 and 8) under the null hypothesis that the
#' sensitivity correlation is zero. Datasets are sampled from this null model
#' and allow for an empirical p-value to be computed that is only significant
#' if the sensitivity correlation is higher than can be expected by chance
#' given the number of samples, correlation and number of miRNAs. p-values
#' are adjusted indepdenently for each parameter
#' combination using Benjamini-Hochberg FDR correction.
#' @export
#'
#' @examples #sponge_compute_p_values()
sponge_compute_p_values <- function(sponge_result,
                                    number.of.samples,
                                    number.of.datasets = 1e5,
                                    simulated.data = NULL){

    #divide gene_gene correlation
    if(max(sponge_result$df) > 7){
        df_breaks <- c(seq(0,7), max(sponge_result$df))
    } else{
        df_breaks <- seq(0, max(sponge_result$df))
    }

    sponge_result <- as.data.table(sponge_result)
    sponge_result <- sponge_result[,
                c("cor_cut", "df_cut") := list(
                    cut(abs(cor), breaks = c(0, seq(0.25, 0.85, 0.1), 1)),
                    cut(df, breaks = df_breaks))]

    levels(sponge_result$cor_cut) <- ks
    levels(sponge_result$df_cut) <- ms

    isplitDT2 <- function(x, ks, ms) {
        ival <- iter(apply(expand.grid(ks, ms), 1, list))
        nextEl <- function() {
            val <- nextElem(ival)
            list(value=x[.(as.character(val[[1]][1]),
                           as.character(val[[1]][2]))], key=val[[1]])
        }
        obj <- list(nextElem=nextEl)
        class(obj) <- c('abstractiter', 'iter')
        obj
    }

    dtcomb <- function(...) {
        rbindlist(list(...))
    }

    setkey(sponge_result, cor_cut, df_cut)

    result <- foreach(dt.m=isplitDT2(sponge_result, ks, ms),
        .combine='dtcomb',
        .multicombine=TRUE,
        .export = c("compute_p_values",
                    "sample_zero_scor_data"),
        .packages = c("gRbase", "MASS", "ppcor", "foreach", "logging", "data.table"),
        .noexport = c("sponge_result")) %dopar% {
            compute_p_values(partition = dt.m$value,
                             cov.matrices = cov.matrices,
                             number.of.datasets = number.of.datasets,
                             number.of.samples = number.of.samples,
                             simulated.data = simulated.data)
        }

    result[p.val == 0, p.val := (1/number.of.datasets)]
    return(result)
}

#' Build null model for p-value computation
#'
#' @param number.of.datasets the number of datesets defining the precision of the p-value
#' @param number.of.samples  the number of samples in the expression data
#'
#' @return a list (for various values of m) of lists (for various values of k)
#' of lists of simulated data sets, drawn from a set of precomputed covariance matrices
#' @import foreach
#' @import data.table
#' @import gRbase
#' @import MASS
#' @import ppcor
#' @export
#'
#' @examples sponge_build_null_model(100, 100)
sponge_build_null_model <- function(number.of.datasets = 1e5,
                                    number.of.samples){
    foreach(cov.matrices.m = cov.matrices[as.character(ms)],
            .final = function(x) setNames(x, as.character(ms)),
            .inorder = TRUE) %:%
        foreach(cov.matrices.k = cov.matrices.m[as.character(ks)],
                .final = function(x) setNames(x, as.character(ks)),
                .inorder = TRUE) %dopar%{
            compute_simulated_data(cov.matrices = cov.matrices.k,
                                   number.of.datasets = number.of.datasets,
                                   number.of.samples = number.of.samples)
        }
}