#' Identify miRNAs for which both genes have miRNA binding sites aka miRNA
#' response elements in the competing endogeneous RNA hypothesis
#'
#' @param geneA The first gene
#' @param geneB The second gene
#' @param mir_interactions A named list of genes, where for each gene all miRNA
#' interacting partners are listed
#'
#' @return A vector with shared RNAs of the two genes.
fn_get_shared_miRNAs <- function(geneA, geneB, mir_interactions){

    source_sign_mirs <- mir_interactions[[geneA]]

    if(!is.null(source_sign_mirs)){
        source_sign_mirs <- as.character(source_sign_mirs$mirna)
    }
    target_sign_mirs <- mir_interactions[[geneB]]

    if(!is.null(target_sign_mirs)){
        target_sign_mirs <- as.character(target_sign_mirs$mirna)
    }
    mir_intersect <- intersect(source_sign_mirs, target_sign_mirs)

    return(mir_intersect)
}

#' Helper function to map MIMAT identifers to miRNA names
#'
#' @param mimats to map
#' @import mirbase.db
#'
#' @return a vector of miRNA names
#'
#' @examples fn_map_mimat_to_mir(c("MIMAT0000086", "MIMAT0000426"))
fn_map_mimats_to_mir <- function(mimats){
    mature_mirnas <- toTable(mirbaseMATURE)
    mature_mirna_names <- mature_mirnas$mature_name
    names(mature_mirna_names) <- mature_mirnas$mature_acc
    mature_mirna_names[mimats]
}

#' Compute all pairwise interactions for a number of genes as indices
#'
#' @param number.of.genes
#' @importFrom gRbase combnPrim
#'
#' @return data frame with one row per unique pairwise combination. To be used
#' as input for the sponge method.
#' @export
#'
#' @examples genes_pairwise_combinations(ncol(gene_expr))
genes_pairwise_combinations <- function(number.of.genes){
    t(combnPrim(number.of.genes, 2))
}

#' Compute competing endogeneous RNA interactions using
#' Sparse Partial correlations ON Gene Expression (SPONGE)
#'
#' @import foreach
#' @import logging
#' @import ppcor
#' @importFrom iterators iter
#' @importFrom iterators icount
#' @importFrom itertools isplitRows
#'
#' @param gene_expr A gene expression matrix
#' @param mir_expr A miRNA expression matrix
#' @param mir_interactions A named list of genes, where for each gene we list
#' all miRNA interaction partners that should be considered.
#' @param log.level The log level, can be one of "info", "debug", "error"
#' @param log.every.n write to the log after every n steps
#' @param selected.genes Operate only on a subset of genes, particularly
#' useful for bootstrapping
#' @param p.adj.method Multiple testing correction method. see ?p.adjust for
#' details
#' @param p.value.threshold Multiple testing p-value cutoff
#' @param parallel.chunks Split into this number of tasks if parallel processing
#' is set up. The number should be high enough to guarantee equal distribution
#' of the work load in parallel execution. However, if the number is too large,
#' e.g. in the worst case one chunk per computation, the overhead causes more
#' computing time than can be saved by parallel execution. Register a parallel
#' backend that is compatible with foreach to use this feature. More information
#' can be found in the documentation of the foreach / doParallel packages.
#'
#' @return A data frame with significant gene-gene competetive endogenous RNA
#' or 'sponge' interactions
#' @export
#'
#' @examples
#' #First we extract miRNA candidates for each of the genes
#' genes_miRNA_candidates <- gene_miRNA_interaction_filter(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' binding_db = NULL)
#'
#' #Second we compute ceRNA interactions for all pairwise combinations of genes
#' #using all miRNAs remaining after filtering through elasticnet.
#' ceRNA_interactions <- sponge(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' mir_interactions = genes_miRNA_candidates)
sponge <- function(gene_expr, mir_expr, mir_interactions,
                  log.level = "INFO",
                  log.every.n = 1e5,
                  selected.genes = NULL,
                  genes.pairwise.combinations = NULL,
                  p.adj.method = "BH",
                  p.value.threshold = 0.05,
                  parallel.chunks = 1e3){
    #names of genes for which we have expr and miRNA data
    genes <- intersect(colnames(gene_expr), names(mir_interactions))

    #now only compute for selected genes
    if(is.null(selected.genes)){
        sel.genes <- genes
    }
    else{
        available.selected.genes <- intersect(selected.genes, genes)
        if(length(available.selected.genes) == 0){
            stop("None of the selected genes is found in the data")
        }
        else if(length(available.selected.genes) < length(selected.genes)){
            warning(paste("Some genes are not found in the data:",
                  paste(setdiff(selected.genes, genes), collapse=","), sep=""))
        }
        sel.genes <- available.selected.genes
    }

    #consider only genes that have miRNA interactions
    sel.genes <- Filter(Negate(is.null), sel.genes)

    #all pairwise combinations of selected genes
    if(is.null(genes.pairwise.combinations)){
        loginfo("Computing all pairwise combinations of genes")
        genes.pairwise.combinations <-
           genes_pairwise_combinations(length(sel.genes))
    }

    loginfo("Beginning SPONGE run...")

    if(is.null(mir_interactions)){
        logwarn("No information on miRNA gene interactions was provided,
                all miRNAs will be considered and runtime will likely explode.")
        mir_intersect <- colnames(mir_expr)
    }
    num_of_samples <- nrow(gene_expr)
    num_of_tasks <- min(nrow(genes.pairwise.combinations), parallel.chunks)

    SPONGE_result <- foreach(gene_combis =
                    itertools::isplitRows(genes.pairwise.combinations,
                                            chunks = parallel.chunks),
                    i = iterators::icount(),
                    .combine=rbind,
                    .packages = c("logging", "ppcor", "foreach", "iterators"),
                    .export = c("fn_get_shared_miRNAs")) %dopar% {
        basicConfig(level = log.level)

        loginfo(paste("SPONGE: worker is processing chunk: ", i, sep=""))

        result <- foreach(gene_combi = iter(gene_combis, by="row"),
                            .combine=rbind) %do% {

            logdebug(paste("Processing source gene", geneA,
                           "and target gene", geneB ))

            geneA <- sel.genes[gene_combi[1]]
            geneB <- sel.genes[gene_combi[2]]

            source_expr <- gene_expr[,geneA]
            target_expr <- gene_expr[,geneB]

            #check if miRNA interaction information is provided, otherwise we
            #consider ALL miRNAs in each comparison
            if(!is.null(mir_interactions)){
                mir_intersect <- fn_get_shared_miRNAs(geneA, geneB,
                                                      mir_interactions)
            }

            #check if there are actually any shared mirnas
            if(length(mir_intersect) == 0){
                logdebug(paste("Source gene", geneA, "and target gene", geneB,
                        "do not share any miRNAs with significant regulation"))
                return(NULL)
            }

            m_expr <- mir_expr[,mir_intersect] #2.c
            pcor <- tryCatch({
                pcor.test(source_expr, target_expr, m_expr)
                #if(length(mir_intersect) == 1){
                    #use linear model here, since glmnet does not
                    #handle single feature models
            #        source_model <- lm(source_expr~m_expr)
            #        target_model <- lm(target_expr~m_expr)

                    #residuals
            #        source_model_residuals <- residuals(source_model) #2f
            #        target_model_residuals <- residuals(target_model) #2f
                #}
                #else{
                #    #lasso
                #    source_model <- cv.glmnet(m_expr, source_expr) #2.d
                #    target_model <- cv.glmnet(m_expr, target_expr) #2.e
                #
                #    #residuals
                #    source_model_residuals <- (predict(source_model, m_expr,
                #    s="lambda.min")  - source_expr)[,1] #2f
                #    target_model_residuals <- (predict(target_model,
                #    m_expr, s="lambda.min")  - target_expr)[,1] #2f
                #}
                #partial correlation
                #cor(source_model_residuals, target_model_residuals)

            }, warning = function(w) {
                logwarn(w)
                return(NULL)
            }, error = function(e) {
                logerror(e)
                return(NULL)
            })

            if(is.null(pcor)) return(NULL)

            #significance
            dof <- pcor$gp #length(mir_intersect)
            dcor <- cor(source_expr, target_expr, use = "complete.obs")
            scor <- dcor - pcor$estimate

            cohens_q <- 0.5 * (log((1+dcor)/(1-dcor)) -
                                   log((1+pcor$estimate)/(1-pcor$estimate)))
            cohens_q_pvalue <- pnorm(cohens_q, lower.tail = F,
                                     sd = sqrt(2 / (num_of_samples - 3)))

            #if(cohens_q_pvalue > p.value.threshold) return(NULL)

            #result
            data.frame(geneA = geneA,
                       geneB = geneB,
                       df = pcor$gp,
                       cor =  dcor,
                       pcor = pcor$estimate,
                       cohens_q = cohens_q,
                       cohens_q_p = cohens_q_pvalue
                      )
        }

        curr_percentage <- round((i / num_of_tasks) * 100, 2)
        loginfo(paste("SPONGE finished chunk: ", i, "->", curr_percentage,
                      "% completed.", sep=""))

        return(result)
    }
    #adjust p-values
    #SPONGE_result$pcor_p.adj <- p.adjust(SPONGE_result$pcor_pval,
    #method = p.adj.method)
    #SPONGE_result$scor_p.adj <- p.adjust(SPONGE_result$scor_pval,
    #method = p.adj.method)
    SPONGE_result$cohens_q_p.adj <- p.adjust(SPONGE_result$cohens_q_p,
                                             method = p.adj.method)

    return(SPONGE_result)
}