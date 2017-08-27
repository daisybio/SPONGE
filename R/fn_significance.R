compute_simulated_data <- function(cov_matrices,
                                   number_of_datasets = 1e5,
                                   number_of_samples){

    #to reach the necessary number of datasets we need to find out how many
    #datasets to construct from each covariance matrix we have
    number_of_datasets_per_matrix <- ceiling(number_of_datasets /
                                                 length(cov_matrices))

    mscor <- sample(
        unlist(sample_zero_mscor_data(
            cov_matrices = cov_matrices,
            number_of_datasets = number_of_datasets_per_matrix,
            number_of_samples = number_of_samples)),
            number_of_datasets)

    test_data_dt <- data.table(mscor)
    setkey(test_data_dt, mscor)

    return(test_data_dt)
}

compute_p_values <- function(partition,
                             cov_matrices = NULL,
                             simulated_data = NULL,
                             number_of_datasets = 1e5,
                             number_of_samples = 100){

    if(nrow(partition) == 1 && is.na(partition$geneA))  return(NULL)
    if(is.null(cov_matrices) & is.null(simulated_data))
      stop("either covariance matrices or simulated data have to be provided.")

    if(!("mscor" %in% colnames(partition)))
        stop("sensitivity correlation missing")

    #check which k and m
    k <- as.character(partition[1,cor_cut])
    m <- as.character(partition[1,df_cut])

    loginfo(paste0("Computing p-values for partition m = ", m, " and k = ", k))

    #simulate data using the appropriate covariance matrices
    if(is.null(simulated_data)){

        cov.matrices.partition <- cov_matrices[[m]][[k]]
        if(is.null(cov.matrices.partition)) test_data_dt <- NULL

        else{
            test_data_dt <- compute_simulated_data(
                cov_matrices = cov.matrices.partition,
                number_of_datasets = number_of_datasets,
                number_of_samples = number_of_samples)
        }
    }
    else{
        logdebug(paste0("Using simulation data provided for partition m = ", m,
                        " and k = ", k))
        test_data_dt <- as.data.table(simulated_data)
    }

    if(is.null(test_data_dt)){
        logdebug(paste0("No covariance matrix found for partition m = ", m,
                        " and k = ", k))
        partition$p.val <- NA
        partition$p.adj <- NA
        partition <- as.data.table(partition)
    }
    else{
        logdebug(paste0("Extracting p-values for partition m = ", m,
                        " and k = ", k,
                        "using simulated data from null model."))
        number.of.datasets.on.right.side <- length(test_data_dt$mscor)

        partition$p.val <- (number.of.datasets.on.right.side -
                                test_data_dt[J(partition$mscor),
                                             .I,
                                             roll = "nearest",
                                             by = .EACHI]$I) /
                                                number.of.datasets.on.right.side
        partition <- as.data.table(partition)
        partition[, p.adj := p.adjust(p.val, method = "BH")]
    }

    return(partition)
}

ks <- seq(0.2, 0.90, 0.1)
ms <- seq(1, 8, 1)

#' Compute p-values for SPONGE interactions
#'
#' @param sponge_result A data frame from a sponge call
#' @param number_of_samples number of samples in the expression matrices
#' @param number_of_datasets number of datasets to simulate.
#' @param simulated_data optional, pre-computed simulated data
#' @param cov_matrices pre-computed covariance matrices to use in case s
#' simulated data is not already provided via simulated_data.
#' affects the mimimum p-value that can be achieved
#' @param log.level The log level of the logging package
#' @import data.table
#' @import gRbase
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
#' @examples null_model <- sponge_build_null_model(100, 100)
#' sponge_compute_p_values(ceRNA_interactions, 100, 100,
#' simulated_data = null_model)
sponge_compute_p_values <- function(sponge_result,
                                    number_of_samples,
                                    number_of_datasets = 1e5,
                                    cov_matrices = cov.matrices,
                                    simulated_data = NULL,
                                    log.level = "INFO"){

    basicConfig(level = log.level)

    loginfo("Computing empirical p-values for SPONGE results.")
    #divide gene_gene correlation
    if(max(sponge_result$df) > 7){
        df_breaks <- c(seq(0,7), max(sponge_result$df))
    } else{
        df_breaks <- seq(0, max(sponge_result$df))
    }

    sponge_result <- as.data.table(sponge_result)
    sponge_result <- sponge_result[,
                                   c("cor_cut", "df_cut") := list(
                                       cut(abs(cor),
                                       breaks = c(0, seq(0.25, 0.85, 0.1), 1)),
                                       cut(df, breaks = df_breaks))]

    levels(sponge_result$cor_cut) <- ks
    levels(sponge_result$df_cut) <- ms

    isplitDT2 <- function(x, ks, ms) {
        ival <- iter(apply(expand.grid(ks, ms), 1, list))
        nextEl <- function() {
            val <- nextElem(ival)

            list(value=x[.(as.character(val[[1]][1]),
                           as.character(val[[1]][2]))],
                 key=val[[1]],
                 sim.data = simulated_data[[as.character(val[[1]][2])]][[
                     as.character(val[[1]][1])]]
            )
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
                                  "sample_zero_mscor_data"),
                      .packages = c("gRbase", "MASS", "ppcor",
                                    "foreach", "logging", "data.table"),
                      .noexport = c("sponge_result")) %do% {
                          partition <- dt.m$value
                          compute_p_values(
                              partition = partition,
                              cov_matrices = cov_matrices,
                              number_of_datasets = number_of_datasets,
                              number_of_samples = number_of_samples,
                              simulated_data = dt.m$sim.data)
                      }

    result[p.val == 0, p.val := (1/number_of_datasets)]
    result[,cor_cut := NULL]
    result[, df_cut := NULL]
    loginfo("Finished computing p-values.")
    return(result)
}

#' Build null model for p-value computation
#'
#' @param number_of_datasets the number of datesets defining the
#' precision of the p-value
#' @param number_of_samples  the number of samples in the expression data
#' @param cov_matrices pre-computed covariance matrices
#' @param log.level The log level of the logging package
#' @return a list (for various values of m) of lists (for various values of k)
#' of lists of simulated data sets, drawn from a set of precomputed
#' covariance matrices
#' @import foreach
#' @import gRbase
#' @importFrom MASS mvrnorm
#' @import ppcor
#' @import logging
#' @export
#'
#' @examples #sponge_build_null_model(100, 100)
sponge_build_null_model <- function(number_of_datasets = 1e5,
                                    number_of_samples,
                                    cov_matrices = cov.matrices,
                                    log.level = "INFO"){

    foreach(cov.matrices.m = cov_matrices[as.character(ms)],
            m = ms,
            .final = function(x) setNames(x, as.character(ms)),
            .inorder = TRUE) %:%
        foreach(cov.matrices.k = cov.matrices.m[as.character(ks)],
                k = ks,
                .final = function(x) setNames(x, as.character(ks)),
                .inorder = TRUE,
                .packages = c("data.table", "gRbase", "MASS",
                              "ppcor", "logging", "foreach")) %dopar%{
                    if(is.null(cov.matrices.k))
                        stop("Covariance matrix missing for simulating data.")
                    basicConfig(level = log.level)

                    loginfo(
                        paste0(
                            "Simulating data for null model of partition m = ",
                            m, " and k = ", k))

                    compute_simulated_data(
                        cov_matrices = cov.matrices.k,
                        number_of_datasets = number_of_datasets,
                        number_of_samples = number_of_samples)
                }
}