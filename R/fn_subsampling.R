#' Sponge subsampling
#' @importFrom foreach foreach
#'
#' @param subsample.n the number of samples to be drawn in each round
#' @param subsample.repeats how often should the subsampling be done?
#' @param subsample.with.replacement logical, should we allow samples to be used
#' repeatedly
#' @param subsample.plot logical, should the results be plotted as box plots
#' @param gene_expr gene expression matrix as defined in sponge
#' @param mir_expr miRNA expression matrix as defined in sponge
#' @param ... parameters passed on to the sponge function
#'
#' @references sponge
#' @return a summary of the results with mean and standard deviations of the
#' correlation and sensitive correlation.
#' @export
#'
#' @examples sponge_subsampling(gene_expr = gene_expr,
#' mir_expr = mir_expr, mir_interactions = mir_interactions,
#' subsample.n = 10, subsample.repeats = 1)
sponge_subsampling <- function(
                    subsample.n = 100,
                    subsample.repeats = 10,
                    subsample.with.replacement = FALSE,
                    subsample.plot = FALSE,
                    gene_expr,
                    mir_expr,
                    ...){

    subsample_results <-
    foreach(sub.n = subsample.n, .combine = rbind) %do%{
        foreach(r = seq_len(subsample.repeats), .combine = rbind) %do% {
            random_draw <- sample.int(nrow(gene_expr), sub.n, replace = subsample.with.replacement)

            sub_gene_expr <- gene_expr[random_draw,]
            sub_mir_expr <- mir_expr[random_draw,]

            if(subsample.with.replacement){
                rownames(sub_gene_expr) <-
                    make.names(rownames(sub_gene_expr), unique = TRUE)

                rownames(sub_mir_expr) <-
                    make.names(rownames(sub_mir_expr), unique = TRUE)
            }

            result <- sponge(gene_expr = sub_gene_expr,
                             mir_expr = sub_mir_expr, ...)
            result$sub.n <- sub.n
            return(result)
        }
    }

    if(subsample.plot){
        if(!requireNamespace("ggplot2", quietly = TRUE)){
            warning("You need to install the package ggplot2 for plotting.")
        }
        else{
            subsample_mscor_plot <- ggplot2::ggplot(subsample_results,
                      ggplot2::aes(x = paste(geneA, geneB, sep = " - "),
                                   y = cor - pcor)) +
                ggplot2::geom_boxplot(fill = "orange") +
                ggplot2::scale_fill_continuous(name = "",
                        labels = c("multiple miRNA sensitivity correlation")) +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x=
                        ggplot2::element_text(angle=90, hjust=1, vjust = 0.5)) +
                ggplot2::ylab("mscor") +
                ggplot2::xlab("ceRNA interaction")

            if(length(subsample.n) > 1){
                subsample_mscor_plot <- subsample_mscor_plot +
                    ggplot2::facet_wrap(~sub.n, ncol = 1)
            }

            print(subsample_mscor_plot)
        }
    }

    return(subsample_results)
}