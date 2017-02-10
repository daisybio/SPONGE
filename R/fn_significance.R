compute_p_values <- function(partition,
                             cov.matrices,
                             number.of.datasets = 1e6,
                             number.of.samples){

    if(!("scor" %in% colnames(partition))) stop("sensitivity correlation missing")

    #check which k and m
    k <- as.character(partition[1,cor_cut])
    m <- as.character(partition[1,df_cut])

    #simulate data using the appropriate covariance matrices
    cov.matrices.partition <- cov.matrices[[m]][[k]]

    #to reach the necessary number of datasets we need to find out how many
    #datasets to construct from each covariance matrix we have
    number.of.datasets.per.matrix <- ceiling(number.of.datasets / length(cov.matrices.partition))

    scor <- unlist(sample_zero_scor_data(cov.matrices = cov.matrices.partition,
                                         number.of.datasets = number.of.datasets.per.matrix,
                                         number.of.samples = number.of.samples)[1:number.of.datasets])

    test_data_dt <- data.table(scor) #data.table(scor[which(scor > 0)])
    setkey(test_data_dt, scor)
    number.of.datasets.on.right.side <- length(test_data_dt$scor)

    #partition_scor_positive <- partition[which(partition$scor >= 0),]

    partition$p.val <- (number.of.datasets.on.right.side -
        test_data_dt[J(partition$scor),
                     .I,
                     roll = "nearest",
                     by = .EACHI]$I) / number.of.datasets.on.right.side

    return(partition)
}

sponge_compute_p_values <- function(sponge_result,
                                    cov.matrices,
                                    number.of.samples,
                                    number.of.datasets = 10,
                                    ms, ks){

    #divide gene_gene correlation
    if(max(sponge_result$df) > 7) df_breaks <- c(seq(0,7), max(sponge_result$df))
    else df_breaks <- seq(0, max(sponge_result$df))

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

    foreach(dt.m=isplitDT2(sponge_result, ks, ms),
        .combine='dtcomb',
        .multicombine=TRUE,
        .export = c("compute_p_values",
                    "sample_zero_scor_data"),
        .packages = c("gRbase", "MASS", "ppcor", "foreach", "logging", "data.table"),
        .noexport = c("sponge_result")) %dopar% {
            compute_p_values(partition = dt.m$value,
                             cov.matrices = cov.matrices,
                             number.of.datasets = number.of.datasets,
                             number.of.samples = number.of.samples)
        }
}