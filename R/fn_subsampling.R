#' Sponge subsampling
#' @importFrom foreach foreach
#' @import ggplot2
#' @import dplyr
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
#' @examples test <- sponge_subsampling(gene_expr = gene_expr,
#' mir_expr = mir_expr, mir_interactions = mir_interactions)



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
        foreach(r = 1:subsample.repeats, .combine = rbind) %do% {
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
    browser()
    if(subsample.plot){
        subsample_scor_plot <- ggplot(subsample_results,
                                      aes(x = paste(geneA, geneB, sep = " - "),
                                          y = cor - pcor)) +
            geom_boxplot(fill = "orange") +
            #geom_boxplot(aes(y = cor), fill = "grey") +
            scale_fill_continuous(name = "",
                    labels = c("multiple miRNA sensitivity correlation")) +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) +
            ylab("mscor") +
            xlab("ceRNA interaction")

        if(length(subsample.n) > 1){
            subsample_scor_plot <- subsample_scor_plot +
                facet_wrap(~sub.n, ncol = 1)
        }

        print(subsample_scor_plot)
    }

    subsample_results %>% dplyr::group_by(geneA, geneB, df) %>%
        dplyr::summarize(cor_mean = mean(cor),
                  cor_sd = sd(cor),
                  pcor_mean = mean(pcor),
                  pcor_sd = sd(pcor))

}