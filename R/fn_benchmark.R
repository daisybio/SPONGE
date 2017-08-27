#' run sponge benchmark where various settings, i.e. with or without
#' regression, single or pooled miRNAs, are compared.
#'
#' @param gene_expr a gene expression matrix
#' @param mir_expr a miRNA expression matrix
#' @param mir_predicted_targets (a list of) mir interaction sources such as
#' targetscan, mircode, etc.
#' @param number_of_genes_to_test a vector of numbers of genes to be tested,
#' e.g. c(250,500)
#' @param number_of_samples number of samples in the null model
#' @param number_of_datasets number of datasets to sample from the null model
#' @param folder where the results should be saved
#'
#' @return a list (regression, no regression) of lists (single miRNA,
#' pooled miRNAs) of benchmark results
#' @export
#'
#' @import logging
#' @import foreach
#'
#' @examples #sponge_run_benchmark(gene_expr = gene_expr, mir_expr = mir_expr,
#' #mir_predicted_targets = targetscan_symbol,
#' #number_of_genes_to_test = c(10), folder = NULL)
sponge_run_benchmark <- function(gene_expr,
                                 mir_expr,
                                 mir_predicted_targets = targetscan_ensg,
                                 number_of_samples = 100,
                                 number_of_datasets = 1e2,
                                 number_of_genes_to_test = c(25),
                                 folder = ""){
    basicConfig(level = "INFO")

    for(num_of_genes in number_of_genes_to_test){
        loginfo(paste("benchmarking with", num_of_genes, "genes"))

        gene_expr_sample <- gene_expr[,sample(colnames(gene_expr),
                                              num_of_genes)]

        gene_miRNA_interaction_results <- foreach(
            elastic.net = c(TRUE, FALSE),
            .final = function(x) setNames(x, c("regression", "no regression")),
            .inorder = TRUE) %do% {
                loginfo(
                    paste(
                        "computing miRNA-gene interactions with elastic.net =",
                        elastic.net))
                miRNA_interactions_time <- system.time(
                    genes_miRNA_candidates <-
                        sponge_gene_miRNA_interaction_filter(
                            gene_expr = gene_expr_sample,
                            mir_expr = mir_expr,
                            elastic.net = elastic.net,
                            mir_predicted_targets = mir_predicted_targets,
                            coefficient.threshold =  -0.05
                        )
                )
                attr(genes_miRNA_candidates, "cputime") <-
                    sum(miRNA_interactions_time[c(1,2,4,5)])
                attr(genes_miRNA_candidates, "elapsedtime") <-
                    miRNA_interactions_time[3]
                return(genes_miRNA_candidates)
            }

        sponge_results <- foreach(
            elastic.net = c("regression", "no regression"),
            .final = function(x) setNames(x, c("regression", "no regression")),
            .inorder = TRUE) %do% {
                foreach(
                    each.miRNA = c(TRUE, FALSE),
                    .final = function(x){
                            setNames(x, c("single miRNA", "pooled miRNAs"))
                        },
                    .inorder = TRUE) %do% {

                        loginfo(paste(
                            "computing miRNA-gene interactions with elastic.net =",
                                      elastic.net, "and considering",
                                      each.miRNA))

                        sponge_time <- system.time(
                            sponge_result <-
                                sponge(
                                    gene_expr = gene_expr_sample,
                                    mir_expr = mir_expr,
                                    mir_interactions =
                                        gene_miRNA_interaction_results[[elastic.net]],
                                    each.miRNA = each.miRNA)
                        )

                        significance_time <- system.time(
                            sponge_result_sign <- sponge_compute_p_values(
                                sponge_result = sponge_result, simulated_data =
                                    NULL, number_of_samples = number_of_samples,
                                number_of_datasets = number_of_datasets) )

                        attr(sponge_result_sign, "cputime_wo_pval") <-
                            sum(sponge_time[c(1,2,4,5)])
                        attr(sponge_result_sign, "elapsedtime_wo_pval") <-
                            sponge_time[3]
                        attr(sponge_result_sign, "cputime") <-
                            sum(sponge_time[c(1,2,4,5)]) +
                            sum(significance_time[c(1,2,4,5)])
                        attr(sponge_result_sign, "elapsedtime") <-
                            sponge_time[3] + significance_time[3]

                        return(sponge_result_sign)
                    }
            }

        if(!is.null(folder)){
            start_date <- date()
            save(sponge_results,
                 gene_miRNA_interaction_results,
                 gene_expr_sample,
                 file = paste(folder,
                              "/benchmark_result_",
                              num_of_genes,
                              "_genes_",
                              str_replace_all(start_date, " ", "_"),
                              ".Rdata",
                              sep = ""))

        }

        return(sponge_results)
    }

}


