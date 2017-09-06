compute_null_model <- function(cov_matrices,
                                   number_of_datasets = 1e5,
                                   number_of_samples){

    #to reach the necessary number of datasets we need to find out how many
    #datasets to construct from each covariance matrix we have
    number_of_datasets_per_matrix <- ceiling(number_of_datasets /
                                                 length(cov_matrices))

    precomputed_cov_matrices <- cov_matrices

    mscor <- sample(
        unlist(sample_zero_mscor_data(
            cov_matrices = precomputed_cov_matrices,
            number_of_datasets = number_of_datasets_per_matrix,
            number_of_samples = number_of_samples)),
            number_of_datasets)

    test_data_dt <- data.table(mscor)
    setkey(test_data_dt, mscor)

    return(test_data_dt)
}

compute_p_values <- function(partition,
                             null_model,
                             number_of_datasets){

    if(!("mscor" %in% colnames(partition)))
        stop("sensitivity correlation missing")

    #check which k and m
    k <- as.character(partition[1,cor_cut])
    m <- as.character(partition[1,df_cut])

    logdebug(paste0("Computing p-values for partition m = ", m, " and k = ", k))

    #simulate data using the appropriate covariance matrices
    logdebug(paste0("Using simulation data provided for partition m = ", m,
                        " and k = ", k))
    test_data_dt <- as.data.table(null_model)

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
        partition[p.val == 0, p.val := (1/number_of_datasets)]
        partition[, p.adj := p.adjust(p.val, method = "BH")]
    }

    return(partition)
}

#given that we can not compute null models for every parameter combination of
#gene-gene correlation k and number of miRNAs m, we need to assign each ceRNA
#interaction in sponge_result to its closest matching null models. we thus
#need a series of ks and ms and assign each ceRNA interaction to its closest
#matching null model.
determine_cutoffs_for_null_model_partitioning <- function(sponge_result,
                                                          ks,
                                                          m_max) {

    if(any(!c("df", "cor") %in% colnames(sponge_result)))
        stop("parameter sponge_result is missing expected columns cor and df")
    if(!is.numeric(sponge_result$df)) stop("column df is not numeric")
    if(!is.numeric(sponge_result$cor)) stop("column df is not numeric")
    if(any(ks <= 0) | any(ks >= 1)) stop("elements in ks outside 0 < k < 1")
    if(m_max < 1) stop("m_max has to be >= 1")

    sponge_result <- as.data.table(sponge_result)
    ms <- seq_len(m_max)

    #set breaks right into the middle of two elements and add 0 and 1 at the
    #boundaries
    cor_breaks <- c(0, ks[-length(ks)] + ((ks[-1] - ks[-length(ks)]) / 2), 1)

    #set partition for m > m_max to m_max and m otherwise
    if(max(sponge_result$df) > (m_max - 1)){
        df_breaks <- c(seq(0,(m_max -1 )), max(sponge_result$df))
    } else{
        df_breaks <- seq(0, max(sponge_result$df))
    }

    sponge_result <- sponge_result[,
                                   c("cor_cut", "df_cut") := list(
                                       cut(abs(cor),
                                           breaks = cor_breaks),
                                       cut(df, breaks = df_breaks))]

    levels(sponge_result$cor_cut) <- ks
    levels(sponge_result$df_cut) <- ms

    setkey(sponge_result, cor_cut, df_cut)

    return(sponge_result)
}

#function for iterating through partitions
isplitDT2 <- function(x, ks, ms, null_model) {
    ival <- iter(apply(expand.grid(ks, ms), 1, list))
    nextEl <- function() {
        val <- nextElem(ival)

        k <- val[[1]][1]
        m <- val[[1]][2]

        sim_data <- null_model[[as.character(m)]][[
            as.character(k)]]

        if(is.null(sim_data))
            stop(paste0("simulation data missing for partition k = ", k,
                        " and m = ", m))

        value <- x[.(as.character(k),
                     as.character(m))]

        if(nrow(value) == 1) if(is.na(value$geneA)) value <- NULL
        list(value = value,
             key = val[[1]],
             sim.data = sim_data
        )
    }
    obj <- list(nextElem=nextEl)
    class(obj) <- c('abstractiter', 'iter')
    obj
}

#function for merging data tables efficiently
dtcomb <- function(...) {
    rbindlist(list(...))
}


#' Compute p-values for SPONGE interactions
#'
#' @param sponge_result A data frame from a sponge call
#' @param null_model optional, pre-computed simulated data
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
#' sponge_compute_p_values(ceRNA_interactions,
#' null_model = null_model)
sponge_compute_p_values <- function(sponge_result,
                                    null_model,
                                    log.level = "ERROR"){

    if(length(null_model) == 0) stop("null model seems to be empty")
    ks <- names(null_model[[1]])
    ms <- names(null_model)

    basicConfig(level = log.level)

    loginfo("Computing empirical p-values for SPONGE results.")

    sponge_result <-
        determine_cutoffs_for_null_model_partitioning(
            sponge_result,
            ks = as.numeric(as.character(ks)),
            m_max = max(as.integer(as.character(ms))))

    number_of_datasets <- nrow(null_model[[1]][[1]])

    result <- foreach(dt.m=isplitDT2(sponge_result, ks, ms, null_model),
                      .combine='dtcomb',
                      .multicombine=TRUE,
                      .export = c("compute_p_values",
                                  "sample_zero_mscor_data"),
                      .packages = c("gRbase", "MASS", "ppcor",
                                    "foreach", "logging", "data.table"),
                      .noexport = c("sponge_result")) %dopar% {
                          partition <- dt.m$value
                          if(is.null(partition)) return(NULL)
                          compute_p_values(
                              partition = partition,
                              null_model = dt.m$sim.data,
                              number_of_datasets = number_of_datasets)
                      }

    result[,cor_cut := NULL]
    result[, df_cut := NULL]
    loginfo("Finished computing p-values.")
    return(as.data.frame(result))
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
#' @examples sponge_build_null_model(100, 100)
sponge_build_null_model <- function(number_of_datasets = 1e5,
                                    number_of_samples,
                                    cov_matrices = precomputed_cov_matrices,
                                    ks = seq(0.2, 0.90, 0.1),
                                    m_max = 8,
                                    log.level = "ERROR"){

    loginfo("Constructing SPONGE null model.")

    if(number_of_datasets < 1) stop("number_of_datasets has to be >= 1")
    if(any(ks <= 0) | any(ks >= 1)) stop("all ks have to be >0 and <1")
    if(m_max < 1) stop("m_max has to be >= 1")

    ms <- seq_len(m_max)
    if((number_of_samples - 2 - m_max) <= 1)
        stop(paste0("sample number to small for m_max = ", m_max))

    null_model <- foreach(cov.matrices.m = cov_matrices[as.character(ms)],
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

                    logdebug(
                        paste0(
                            "Simulating data for null model of partition m = ",
                            m, " and k = ", k))

                    compute_null_model(
                        cov_matrices = cov.matrices.k,
                        number_of_datasets = number_of_datasets,
                        number_of_samples = number_of_samples)
                              }
    loginfo("Finished constructing SPONGE null model.")
    return(null_model)
}