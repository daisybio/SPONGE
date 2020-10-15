#iterate over chunks of rows for efficient parallel computation
split_cols <- function(x, ...)
{
  it <- idiv(ncol(x), ...)
  i <- 1L
  nextEl <- function() {
    n <- as.integer(nextElem(it))
    j <- i
    i <<- i + n
    x[, seq(j, length = n) , drop = FALSE]
  }
  object <- list(nextElem = nextEl)
  class(object) <- c("abstractiter", "iter")
  object
}


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
#' @import doRNG
#' @importFrom iterators iter
#'
#' @param gene_expr A gene expression matrix with samples in rows and featurs
#' in columns. Alternatively an object of class ExpressionSet.
#' @param mir_expr A miRNA expression matrix with samples in rows and features
#' in columns. Alternatively an object of class ExpressionSet.
#' @param mir_predicted_targets A data frame with miRNA in cols and genes in rows.
#' A 0 indicates the miRNA is not predicted to target the gene, >0 otherwise.
#' If this parameter is NULL all miRNA-gene interactions are tested
#' @param log.file Log file to write to
#' @param log.level One of 'warn', 'error', 'info'
#' @param var.threshold Only consider genes and miRNA with
#' variance > var.threshold. If this parameter is NULL no variance filtering
#' is performed.
#' @param random_seed A random seed to be used for reproducible results
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
#' @param parallel.chunks Split into this number of tasks if parallel processing
#' is set up. The number should be high enough to guarantee equal distribution
#' of the work load in parallel execution. However, if the number is too large,
#' e.g. in the worst case one chunk per computation, the overhead causes more
#' computing time than can be saved by parallel execution. Register a parallel
#' backend that is compatible with foreach to use this feature. More information
#' can be found in the documentation of the foreach / doParallel packages.
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
                                                 log.level = "ERROR",
                                                 log.file = NULL,
                                                 var.threshold = NULL,
                                                 F.test = FALSE,
                                                 F.test.p.adj.threshold = 0.05,
                                                 coefficient.threshold = -0.05,
                                                 coefficient.direction = "<",
                                                 select.non.targets = FALSE,
                                                 random_seed = NULL,
                                                 parallel.chunks = 100){
    if(!is.null(log.file))
        addHandler(writeToFile, file=log.file, level=log.level)
    else{
        basicConfig(level = log.level)
    }
    with_target_info <- !is.null(mir_predicted_targets)
    foreach_packages <- c("logging", "glmnet")

    if(is(gene_expr, "big.matrix.descriptor") && requireNamespace("bigmemory"))
    {
        loginfo("Detected gene expression big matrix descriptor")
        gene_expr <- check_and_convert_expression_data(gene_expr)
        gene_expr <- bigmemory::as.matrix(gene_expr)
    }
    else{
        gene_expr <- check_and_convert_expression_data(gene_expr)
    }

    if(is(mir_expr, "big.matrix.descriptor") && requireNamespace("bigmemory"))
    {
        loginfo("Detected miRNA expression big matrix descriptor")
        mir_expr <- check_and_convert_expression_data(mir_expr)
        mir_expr <- bigmemory::as.matrix(mir_expr)
    }
    else{
        mir_expr <- check_and_convert_expression_data(mir_expr)
    }

    if(select.non.targets)
        logwarn("Selecting only miRNA targets not predicted as targets")

    if(!is.matrix(gene_expr) | any(dim(gene_expr) < 2)) stop("gene_expr matrix not properly formatted")
    if(!is.matrix(gene_expr) | any(dim(mir_expr) < 2)) stop("mir_expr matrix not properly formatted")
    if(nrow(mir_expr) != nrow(gene_expr))
        stop("mir_expr and gene_expr matrix differ in row numbers")

    #remove columns with little variance
    if(!is.null(var.threshold)){
        loginfo("Removing genes and miRNAs below variance threshold")
        gene_expr <- gene_expr[,which(apply(gene_expr, 2, var) > var.threshold)]
        mir_expr <- mir_expr[,which(apply(mir_expr, 2, var) > var.threshold)]
    } else{
        loginfo("Removing genes and miRNAs with zero variance")
        gene_expr <- gene_expr[,which(apply(gene_expr, 2, var) != 0)]
        mir_expr <- mir_expr[,which(apply(mir_expr, 2, var) != 0)]
    }

    if(!is.matrix(gene_expr) | any(dim(gene_expr) < 2))
        stop("variance threshold too strict in gene_expr")
    if(!is.matrix(mir_expr) | any(dim(mir_expr) < 2))
        stop("variance threshold too strict in mir_expr")

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


        mir_predicted_targets <- matrix(
          nrow = length(all_genes),
          ncol = length(all_mirs),
          dimnames = list(all_genes, all_mirs))

        mir_predicted_targets[,] <- 0

        #fill matrix
        for(mir_db in list_of_mir_predicted_targets){
            mir_db_genes <- which(all_genes %in% rownames(mir_db))
            mir_db_mirs <- which(all_mirs %in% colnames(mir_db))
            mir_predicted_targets[mir_db_genes, mir_db_mirs] <-
                mir_predicted_targets[mir_db_genes, mir_db_mirs] +
                mir_db[all_genes[mir_db_genes], all_mirs[mir_db_mirs]]
        }
    }
    else{
        all_mirs <- intersect(colnames(mir_predicted_targets), colnames(mir_expr))
        all_genes <- intersect(rownames(mir_predicted_targets), colnames(gene_expr))

        mir_predicted_targets <- mir_predicted_targets[all_genes, all_mirs]
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

    num_of_tasks <- min(max(1, ceiling(ncol(gene_expr) / 10)), parallel.chunks)

    loginfo("Computing gene / miRNA regression models...")

    #loop over all genes and compute regression models to identify important miRNAs
    final_result <- foreach(g_expr_batch = split_cols(gene_expr, chunks = num_of_tasks),
            chunk = seq_len(num_of_tasks),
            .packages = foreach_packages,
            .export = c("fn_get_model_coef", "fn_elasticnet", "fn_gene_miRNA_F_test", "fn_get_rss"),
            .inorder = TRUE,
            .options.RNG = random_seed
            ) %dorng% {

              #setup logging
                if(!is.null(log.file))
                    addHandler(writeToFile, file=log.file, level=log.level)
                else{
                    basicConfig(level = log.level)
                }
              loginfo(paste("Computing gene / miRNA regression models: chunk ", chunk, " of ", num_of_tasks, ".", sep=""))

              batch_result <- foreach(g_expr = iterators::iter(g_expr_batch, by = "col"),
                                      gene = colnames(g_expr_batch),
                                      .final = function(x) setNames(x, colnames(g_expr_batch))) %do%{

                gene_idx <- which(all_genes == gene)

                #check if we have a target database
                if(with_target_info){
                    mimats_matched <- all_mirs[which(mir_predicted_targets[gene_idx,] > 0)]

                    if(length(mimats_matched) == 0){
                        return(NULL)
                    }

                    if(select.non.targets){
                        non_targets <- all_mirs[which(mir_predicted_targets[gene_idx,] == 0)]
                        mimats_matched <- sample(non_targets,
                                                 min(length(mimats_matched), length(non_targets)))
                    }
                    m_expr <- mir_expr[,which(all_mirs %in% mimats_matched)]
                }
                else{
                    mimats_matched <- all_mirs
                    m_expr <- as.matrix(mir_expr)
                }

                if(!elastic.net){
                    return(data.frame(mirna = mimats_matched))
                }

                #learn a regression model to figure out which miRNAs regulate this gene in
                #the given dataset

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
              return(batch_result)
            }
    loginfo("FINISHED")
    return(unlist(final_result, recursive = FALSE))
}
