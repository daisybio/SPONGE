#' Determine miRNA-gene interactions to be considered in SPONGE
#' @description The purpose of this method is to limit the number of miRNA-gene
#' interactions we need to consider in SPONGE. There are 3 filtering steps:
#' 1. variance filter (optional). Only considre genes and miRNAs with variance >
#' var.threshold.
#' 2. miRNA target database filter (optional). Use a miRNA target database
#' provided by the user to filter for those miRNA gene interactions for which
#' evidence exists. This can either be predicted target interactions or
#' experimentally validated ones.
#' 3. For each remaining interaction of a gene and its regulating miRNAs
#' use elastic net regression to achieve
#' a) Feature selection: We only retain miRNAs that influence gene expression
#' b) Effect strength: The sign of the coefficients allows
#' us to filter for miRNAs that down-regulate gene expression. Moreover, we can
#' use the coefficients to rank the miRNAs by their relative effect strength.
#' We strongly recommend setting up a parallel backend compatible with the
#' foreach package. See example and the documentation of the
#' foreach and doParallel packages.
#' @import foreach
#' @import glmnet
#' @import bigmemory
#'
#' @param gene_expr A matrix of gene expression
#' @param mir_expr A matrix of miRNA expression
#' @param mir_predicted_targets A data frame with miRNA in cols and genes in rows.
#' A 0 indicates the miRNA is not predicted to target the gene, >0 otherwise.
#' If this parameter is NULL all miRNA-gene interactions are tested
#' @param log.every.n Log only ever n iterations to limit output
#' @param log.level One of 'warn', 'error', 'info'
#' @param var.threshold Only consider genes and miRNA with
#' variance > var.threshold. If this parameter is NULL no variance filtering
#' is performed.
#' @param F.test If true, an F-test is performed on each model parameter to
#' assess its importance for the model based on the RSS of the full model vs
#' the RSS of the nested model without the miRNA in question. This is time
#' consuming and has the potential disadvantage that correlated miRNAs are
#' removed even though they might play a role in ceRNA interactions. Use at your
#' own risk.
#' @param F.test.p.adj.threshold If F.test is TRUE, threshold to use for
#' miRNAs to be included.
#' @param coefficient.direction If "<", coefficient has to be lower than
#' coefficient.threshold, if ">", coefficient has to be larger than threshold.
#' If NULL, the absolute value of the coefficient has to be larger than the
#' threshold.
#' @param coefficient.threshold threshold to cross for a regression coefficient
#' to be called significant. depends on the parameter coefficient.direction.
#' @param select.non.targets For testing effect of miRNA target information.
#' If TRUE, the method determines as usual which miRNAs are potentially
#' targeting a gene. However, these are then replaced by a random sample of
#' non-targeting miRNAs (without seeds) of the same size. Useful for testing
#' if observed effects are caused by miRNA regulation.
#' @param elastic.net Whether to apply elastic net regression filtering or not.
#' @return A list of genes, where for each gene, the regulating miRNA are
#' included as a data frame. For F.test = TRUE this is a data frame with fstat
#' and p-value for each miRNA. Else it is a data frame with the model
#' coefficients.
#' @export
#'
#' @seealso sponge
#' @examples
#' #library(doParallel)
#' #cl <- makePSOCKcluster(2)
#' #registerDoParallel(cl)
#' genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' mir_predicted_targets = targetscan_symbol)
#' #stopCluster(cl)
#'
#' #If we also perform an F-test, only few of the above miRNAs remain
#' genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' mir_predicted_targets = targetscan_symbol,
#' F.test = TRUE,
#' F.test.p.adj.threshold = 0.05)
#'
sponge_gene_miRNA_interaction_filter <- function(gene_expr, mir_expr,
                                                 mir_predicted_targets,
                                                 elastic.net = TRUE,
                                                 log.every.n = 100,
                                                 log.level = "OFF",
                                                 var.threshold = NULL,
                                                 F.test = FALSE,
                                                 F.test.p.adj.threshold = 0.05,
                                                 coefficient.threshold = -0.05,
                                                 coefficient.direction = "<",
                                                 select.non.targets = FALSE){
    basicConfig(level = log.level)
    with_target_info <- !is.null(mir_predicted_targets)

    #remove columns with little variance
    if(!is.null(var.threshold)){
        loginfo("Removing genes and miRNAs below variance threshold")
        gene_expr <- gene_expr[,which(colVars(gene_expr) > var.threshold)]
        mir_expr <- mir_expr[,which(colVars(mir_expr) > var.threshold)]
    } else{
        loginfo("Removing genes and miRNAs with zero variance")
        gene_expr <- gene_expr[,which(apply(gene_expr, 2, var) != 0)]
        mir_expr <- mir_expr[,which(apply(mir_expr, 2, var) != 0)]
    }

    #merge mirna target annotation
    loginfo("merging miRNA target database annotations")
    if(!with_target_info){
        all_mirs <- colnames(mir_expr)
        all_genes <- colnames(gene_expr)
    }
    else if(is.list(mir_predicted_targets)){
        list_of_mir_predicted_targets <- mir_predicted_targets

        all_mirs <- foreach(mir_db = mir_predicted_targets,
                            .combine = union) %do% {
                                colnames(mir_db)
                            }
        all_mirs <- intersect(all_mirs, colnames(mir_expr))

        all_genes <- foreach(mir_db = mir_predicted_targets,
                             .combine = union) %do% {
                                 rownames(mir_db)
                             }
        all_genes <- intersect(all_genes, colnames(gene_expr))

        mir_predicted_targets <- big.matrix(
            nrow = length(all_genes),
            ncol = length(all_mirs),
            dimnames = list(all_genes, all_mirs))

        mir_predicted_targets[,] <- 0

        #fill big matrix
        for(mir_db in list_of_mir_predicted_targets){
            mir_db_genes <- intersect(all_genes, rownames(mir_db))
            mir_db_mirs <- intersect(all_mirs, colnames(mir_db))
            mir_predicted_targets[mir_db_genes, mir_db_mirs] <-
                mir_predicted_targets[mir_db_genes, mir_db_mirs] +
                mir_db[mir_db_genes, mir_db_mirs]
        }

        mir_predicted_targets_description <- describe(mir_predicted_targets)
    }
    else{
        all_mirs <- intersect(colnames(mir_predicted_targets), colnames(mir_expr))
        all_genes <- intersect(rownames(mir_predicted_targets), colnames(gene_expr))

        mir_predicted_targets <- mir_predicted_targets[all_genes, all_mirs]
        mir_predicted_targets <- as.big.matrix(mir_predicted_targets)

        mir_predicted_targets_description <- describe(mir_predicted_targets)
    }

    omitted_mirnas <- setdiff(colnames(mir_expr), all_mirs)
    omitted_genes <- setdiff(colnames(gene_expr), all_genes)

    if(length(omitted_mirnas) > 0){
        logwarn(paste0(length(omitted_mirnas), " miRNAs were omitted because we do not have miRNA target interaction data for them"))
        logdebug(paste0("miRNAs ",paste(omitted_mirnas, collapse = "/"), " were omitted because we do not have miRNA target interaction data for them"))
    }

    if(length(omitted_genes) > 0){
        logwarn(paste0(length(omitted_genes), " genes were omitted because we do not have miRNA target interaction data for them"))
        logdebug(paste0("genes ",paste(omitted_genes, collapse = "/"), " were omitted because we do not have miRNA target interaction data for them"))
    }
    gene_expr <- gene_expr[,all_genes]
    mir_expr <- mir_expr[,all_mirs]
    num_of_genes <- length(all_genes)

    #convert to big.matrix for shared memory operations
    gene_expr_bm <- as.big.matrix(gene_expr)
    gene_expr_description <- describe(gene_expr_bm)

    mir_expr_bm <- as.big.matrix(mir_expr)
    mir_expr_description <- describe(mir_expr_bm)

    #loop over all genes and compute regression models to identify important miRNAs
    foreach(gene_idx = seq_len(num_of_genes),
            .final = function(x) setNames(x, all_genes),
            .packages = c("logging", "glmnet", "dplyr", "bigmemory"),
            .export = c("fn_get_model_coef", "fn_elasticnet", "fn_gene_miRNA_F_test", "fn_get_rss"),
            .inorder = TRUE) %dopar% {

                #attach shared data
                attached_gene_expr <- attach.big.matrix(gene_expr_description)
                attached_mir_expr <- attach.big.matrix(mir_expr_description)

                #get gene name
                gene <- all_genes[gene_idx]

                #setup logging
                basicConfig(level = log.level)

                if(gene_idx %% log.every.n == 0){
                    curr_percentage <- round((gene_idx / num_of_genes) * 100, 2)
                    loginfo(paste("Computing gene / miRNA regression models: ", curr_percentage, "% completed.", sep=""))
                }

                logdebug(paste("Processing gene", gene))

                #expression values of this gene
                g_expr <- attached_gene_expr[,gene_idx]

                #check if we have a target database
                if(with_target_info){
                    attached_mir_predicted_targets <- attach.big.matrix(mir_predicted_targets_description)
                    mimats_matched <- all_mirs[which(attached_mir_predicted_targets[gene_idx,] > 0)]

                    if(length(mimats_matched) == 0){
                        logdebug("None of the target mirnas are found in expression data. Returning null for this gene.")
                        return(NULL)
                    }

                    if(select.non.targets){
                        non_targets <- all_mirs[which(attached_mir_predicted_targets[gene_idx,] == 0)]
                        mimats_matched <- sample(non_targets,
                                                 min(length(mimats_matched), length(non_targets)))
                    }
                    m_expr <- attached_mir_expr[,which(all_mirs %in% mimats_matched)]
                }
                else{
                    mimats_matched <- all_mirs
                    m_expr <- as.matrix(attached_mir_expr)
                }

                if(!elastic.net){
                    return(data.frame(mirna = mimats_matched))
                }

                #learn a regression model to figure out which miRNAs regulate this gene in
                #the given dataset
                logdebug(paste("Learning regression model for gene", gene))

                #if we have only one miRNA we can't use glmnet.
                #We use lm instead and check if this miRNA is sigificant
                if(length(mimats_matched) == 1) {
                    tryCatch({
                        if(F.test){
                            fstat <- as.numeric(summary(lm(g_expr ~ m_expr))$fstatistic[1])
                            pval <- pf(fstat, 1, length(m_expr), lower.tail=FALSE)
                            lm_result <- data.frame(mirna = mimats_matched[1],
                                                    fstats = fstat,
                                                    pval=pval,
                                                    p.adj = pval)
                        }else{
                            lm_result <- data.frame(mirna = mimats_matched[1],
                                                    coefficient = coef(lm(g_expr ~ m_expr))[-1])

                            if(is.null(coefficient.direction) && !is.null(coefficient.threshold))
                                outside.threshold <- which(abs(lm_result$coefficient) > coefficient.threshold)
                            else if(coefficient.direction == "<")
                                outside.threshold <- which(lm_result$coefficient < coefficient.threshold)
                            else if(coefficient.direction == ">")
                                outside.threshold <- which(lm_result$coefficient > coefficient.threshold)

                            if(length(outside.threshold) == 0) return(NULL)
                            else return(lm_result[outside.threshold,])
                        }
                    }, warning = function(w) {
                        logdebug(w)
                        return(NULL)
                    }, error = function(e) {
                        logerror(e)
                        return(NULL)
                    })
                    return(lm_result)
                }

                #if we have more than one miRNA we can use elasticnet to
                #find out which are the essential features
                model <- tryCatch({
                    fn_elasticnet(m_expr, g_expr) #elasticnet trying different alphas
                }, warning = function(w) {
                    logdebug(w)
                    return(NULL)
                }, error = function(e) {
                    logerror(e)
                    return(NULL)
                })
                if(is.null(model)) return(NULL)

                #extract model coefficients
                if(!F.test){
                    result <- fn_get_model_coef(model)

                    if(is.null(coefficient.direction) && !is.null(coefficient.threshold))
                        outside.threshold <- which(abs(result$coefficient) > coefficient.threshold)
                    else if(coefficient.direction == "<")
                        outside.threshold <- which(result$coefficient < coefficient.threshold)
                    else if(coefficient.direction == ">")
                        outside.threshold <- which(result$coefficient > coefficient.threshold)

                    if(length(outside.threshold) == 0) return(NULL)
                    else return(result[outside.threshold,])
                }
                #we use the F test to assess the significance of each feature
                else if(F.test){
                    result <- tryCatch({
                        fn_gene_miRNA_F_test(g_expr, m_expr, model,
                                             F.test.p.adj.threshold)
                    }, warning = function(w) {
                        logdebug(w)
                        return(NULL)
                    }, error = function(e) {
                        logerror(e)
                        return(NULL)
                    })

                    return(result)
                }
            }
}
