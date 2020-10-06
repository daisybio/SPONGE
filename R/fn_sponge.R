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

#' Compute all pairwise interactions for a number of genes as indices
#'
#' @param number.of.genes Number of genes for which all pairwise interactions
#' are needed
#' @importFrom gRbase combnPrim
#'
#' @return data frame with one row per unique pairwise combination. To be used
#' as input for the sponge method.
#'
genes_pairwise_combinations <- function(number.of.genes){
    t(combnPrim(number.of.genes, 2))
}

#' Compute competing endogeneous RNA interactions using
#' Sparse Partial correlations ON Gene Expression (SPONGE)
#'
#' @import foreach
#' @import logging
#' @import doRNG
#' @importFrom ppcor pcor pcor.test
#' @importFrom iterators iter
#' @importFrom iterators icount
#' @importFrom data.table data.table
#' @importFrom gRbase combnPrim
#' @importFrom stats coef cor cov2cor df lm p.adjust pf predict reorder rnorm runif sd setNames xtabs var
#'
#' @param gene_expr A gene expression matrix with samples in rows and featurs
#' in columns. Alternatively an object of class ExpressionSet.
#' @param mir_expr A miRNA expression matrix with samples in rows and features
#' in columns. Alternatively an object of class ExpressionSet.
#' @param mir_interactions A named list of genes, where for each gene we list
#' all miRNA interaction partners that should be considered.
#' @param log.level The log level, can be one of "info", "debug", "error"
#' @param log.every.n write to the log after every n steps
#' @param log.file write log to a file, particularly useful for paralleliyzation
#' @param selected.genes Operate only on a subset of genes, particularly
#' useful for bootstrapping
#' @param gene.combinations A data frame of combinations of genes to be tested.
#' Gene names are taken from the first two columns and have to match the names
#' used for gene_expr
#' @param random_seed A random seed to be used for reproducible results
#' @param result_as_dt whether to return results as data table or data frame
#' @param parallel.chunks Split into this number of tasks if parallel processing
#' is set up. The number should be high enough to guarantee equal distribution
#' of the work load in parallel execution. However, if the number is too large,
#' e.g. in the worst case one chunk per computation, the overhead causes more
#' computing time than can be saved by parallel execution. Register a parallel
#' backend that is compatible with foreach to use this feature. More information
#' can be found in the documentation of the foreach / doParallel packages.
#' @param each.miRNA Whether to consider individual miRNAs or pooling
#' them.
#' @param min.cor Consider only gene pairs with a minimum correlation specified
#' here.
#'
#' @return A data frame with significant gene-gene competetive endogenous RNA
#' or 'sponge' interactions
#' @export
#'
#' @examples
#' #First, extract miRNA candidates for each of the genes
#' #using sponge_gene_miRNA_interaction_filter. Here we use a prepared
#' #dataset mir_interactions.
#'
#' #Second we compute ceRNA interactions for all pairwise combinations of genes
#' #using all miRNAs remaining after filtering through elasticnet.
#' ceRNA_interactions <- sponge(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' mir_interactions = mir_interactions)
sponge <- function(gene_expr,
                   mir_expr,
                   mir_interactions = NULL,
                   log.level = "ERROR",
                   log.every.n = 1e5,
                   log.file = NULL,
                   selected.genes = NULL,
                   gene.combinations = NULL,
                   each.miRNA = FALSE,
                   min.cor = 0.1,
                   parallel.chunks = 1e3,
                   random_seed = NULL,
                   result_as_dt = FALSE){

    if(!is.null(log.file))
        addHandler(writeToFile, file=log.file, level=log.level)
    else{
        basicConfig(level = log.level)
    }

    #handle bigmemory objects and check expression matrices
    foreach_packages <- c("logging", "ppcor", "foreach",
                          "iterators", "data.table")

    if(is(gene_expr, "big.matrix.descriptor") && requireNamespace("bigmemory"))
    {
        loginfo("Detected gene expression big matrix descriptor")
        gene_expr_big_memory <- TRUE
        gene_expr_description <- gene_expr
        gene_expr <- check_and_convert_expression_data(gene_expr)
        foreach_packages <- union(foreach_packages, "bigmemory")
    }
    else{
        gene_expr_big_memory <- FALSE
        gene_expr <- check_and_convert_expression_data(gene_expr)
    }

    if(is(mir_expr, "big.matrix.descriptor") && requireNamespace("bigmemory"))
    {
        loginfo("Detected miRNA expression big matrix descriptor")
        mir_expr_big_memory <- TRUE
        mir_expr_description <- mir_expr
        mir_expr <- check_and_convert_expression_data(mir_expr)
        foreach_packages <- union(foreach_packages, "bigmemory")
    }
    else{
        mir_expr_big_memory <- FALSE
        mir_expr <- check_and_convert_expression_data(mir_expr)
    }

    if(is.null(mir_interactions)){
        logwarn("No information on miRNA gene interactions was provided,
                all miRNAs will be considered and runtime will likely explode.")
        genes <- colnames(gene_expr)
    }
    else{
        #filter out genes without miR interactions
        mir_interactions <- Filter(Negate(is.null), mir_interactions)

        #names of genes for which we have expr and miRNA data
        genes <- intersect(colnames(gene_expr), names(mir_interactions))
    }

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

    #use indices for faster access
    genes.as.indices <- FALSE

    #make sure sel.genes order is used
    gene_expr <- gene_expr[,sel.genes]

    #for getting indices of things
    all_mirs <- colnames(mir_expr)
    all_genes <- colnames(gene_expr)

    #all pairwise combinations of selected genes
    if(is.null(gene.combinations)){
        loginfo("Computing all pairwise combinations of genes")
        #consider only genes that have miRNA interactions

        gene.combinations <-
           genes_pairwise_combinations(length(sel.genes))

        gene.combinations <- data.frame(gene.combinations,
                                   all_genes[gene.combinations[,1]],
                                   all_genes[gene.combinations[,2]],
                                   stringsAsFactors = FALSE)

    } else{
      if(ncol(gene.combinations) > 2)
        stop("gene.combinations is expected to have two columns with gene identifiers")

      colnames(gene.combinations) <- c("geneA", "geneB")

      gene.combinations$geneA <- as.character(gene.combinations$geneA)
      gene.combinations$geneB <- as.character(gene.combinations$geneB)

      valid_genes_A <- which(gene.combinations$geneA %in% all_genes)
      valid_genes_B <- which(gene.combinations$geneB %in% all_genes)
      valid_genes <- intersect(valid_genes_A, valid_genes_B)

      if(length(valid_genes) > 0){
        gene.combinations <- gene.combinations[valid_genes,]

        gene.combinations <- data.frame(
          which(all_genes %in% gene.combinations$geneA),
          which(all_genes %in% gene.combinations$geneB),
          gene.combinations,
          stringsAsFactors = FALSE)
      } else stop("No valid gene combinations selected")
    }
    colnames(gene.combinations) <- c("geneA_idx", "geneB_idx", "geneA", "geneB")

    loginfo("Beginning SPONGE run...")

    if(!gene_expr_big_memory) gene_expr_description <- gene_expr
    if(!mir_expr_big_memory) mir_expr_description <- mir_expr

    rm(gene_expr)
    rm(mir_expr)

    num_of_samples <- nrow(gene_expr)
    num_of_tasks <- min(max(1, ceiling(nrow(gene.combinations) / 1000)),
                        parallel.chunks)

    SPONGE_result <- foreach(
        gene_combis =
            split_rows(gene.combinations,
                       chunks = num_of_tasks),
        i = iterators::icount(),
        .combine = function(...)
            rbindlist(list(...)),
        .multicombine = TRUE,
        .packages = foreach_packages,
        .export = c("fn_get_shared_miRNAs", "processChunk", "compute_pcor"),
        .options.RNG = random_seed
    ) %dorng% {
        if (!is.null(log.file))
            addHandler(writeToFile, file = log.file, level = log.level)
        else{
            basicConfig(level = log.level)
        }

        loginfo(paste("SPONGE: worker is processing chunk: ", i, sep = ""))

        #attach bigmemory objects if necessary. avoid using names of actual
        #matrix objects because they would then be exported to the workers
        if (gene_expr_big_memory)
            attached_gene_expr <- attach.big.matrix(gene_expr_description)
        else
            attached_gene_expr <- gene_expr_description

        if (mir_expr_big_memory)
            attached_mir_expr <- attach.big.matrix(mir_expr_description)
        else
            attached_mir_expr <- mir_expr_description

        #if(require(pryr)) logdebug(paste("current memory used by worker:", pryr::mem_used()))

        result <-
            processChunk(
                gene_combis,
                attached_gene_expr,
                attached_mir_expr,
                mir_interactions,
                all_mirs,
                each.miRNA,
                min.cor
            )

        loginfo(paste("SPONGE finished chunk:", i, "of", num_of_tasks))
        if(is.null(result)) return(list())
        else return(result)
    }
    loginfo("SPONGE completed successfully. Returning results.")

    if(result_as_dt) return(SPONGE_result)
    else return(as.data.frame(SPONGE_result))
}

#internal function
processChunk <- function(gene_combis, attached_gene_expr, attached_mir_expr, mir_interactions,
                         all_mirs, each.miRNA, min.cor){
    if(is.null(mir_interactions))
        mir_intersect <- all_mirs

    foreach(geneA_idx = gene_combis$geneA_idx,
            geneB_idx = gene_combis$geneB_idx,
            geneA = gene_combis$geneA,
            geneB = gene_combis$geneB,
            .export = c("compute_pcor", "fn_get_shared_miRNAs"),
            .combine=function(...) rbindlist(list(...))) %do% {

              source_expr <- attached_gene_expr[,geneA_idx]
              target_expr <- attached_gene_expr[,geneB_idx]

              #check correlation
              dcor <- cor(source_expr, target_expr)

              if(is.na(dcor)) return(NULL)

              if(!is.null(min.cor)){
                if(dcor < min.cor)
                  return(NULL)
              }

              #check if miRNA interaction information is provided, otherwise we
              #consider ALL miRNAs in each comparison
              if(!is.null(mir_interactions)){
                mir_intersect <- fn_get_shared_miRNAs(geneA, geneB,
                                                      mir_interactions)

                #check if shared miRNAs are in expression matrix
                if(length(setdiff(mir_intersect, all_mirs)) > 0){
                  mir_intersect <- intersect(mir_intersect, all_mirs)
                }

                #check if there are actually any shared mirnas
                if(length(mir_intersect) == 0){
                  return(NULL)
                }
              }

              if(each.miRNA){
                result <- foreach(mirna = mir_intersect,
                                  .export = c("compute_pcor"),
                                  .combine = function(...) rbindlist(list(...)),
                                  .inorder = TRUE) %do%{
                            m_expr <- attached_mir_expr[,which(all_mirs == mirna)]
                            compute_pcor(source_expr, target_expr, m_expr,
                                         geneA, geneB, dcor)
                          }
                result$miRNA <- mir_intersect
                return(result)
              }
              else{
                m_expr <- attached_mir_expr[,which(all_mirs %in% mir_intersect)]
                compute_pcor(source_expr, target_expr, m_expr,
                             geneA, geneB, dcor)
              }
            }
}

#iterate over chunks of rows for efficient parallel computation
split_rows <- function(x, ...)
{
    it <- idiv(nrow(x), ...)
    i <- 1L
    nextEl <- function() {
        n <- as.integer(nextElem(it))
        j <- i
        i <<- i + n
        x[seq(j, length = n), , drop = FALSE]
    }
    object <- list(nextElem = nextEl)
    class(object) <- c("abstractiter", "iter")
    object
}

#compute partial correlation
compute_pcor <- function(source_expr, target_expr, m_expr,
                         geneA, geneB, dcor){

    pcor <- tryCatch({
        pcor.test(source_expr, target_expr, m_expr)
    }, warning = function(w) {
        logdebug(w)
        suppressWarnings(pcor.test(source_expr, target_expr, m_expr))
    }, error = function(e) {
        logerror(e)
        return(NULL)
    })

    if(is.null(pcor)) return(NULL)

    list(geneA = geneA,
         geneB = geneB,
         df = pcor$gp,
         cor =  dcor,
         pcor = pcor$estimate,
         mscor = dcor - pcor$estimate
    )
}
