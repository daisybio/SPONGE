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
#' @importFrom ppcor pcor
#' @importFrom iterators iter
#' @importFrom iterators icount
#' @importFrom data.table data.table
#' @importFrom bigmemory describe
#' @importFrom bigmemory as.big.matrix
#' @importFrom bigmemory attach.big.matrix
#' @importFrom gRbase combnPrim
#' @importFrom stats coef cor cov2cor df lm p.adjust pf predict reorder rnorm
#' runif sd setNames xtabs var
#'
#' @param gene_expr A gene expression matrix
#' @param mir_expr A miRNA expression matrix
#' @param mir_interactions A named list of genes, where for each gene we list
#' all miRNA interaction partners that should be considered.
#' @param log.level The log level, can be one of "info", "debug", "error"
#' @param log.every.n write to the log after every n steps
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
                   selected.genes = NULL,
                   gene.combinations = NULL,
                   each.miRNA = FALSE,
                   min.cor = 0.1,
                   parallel.chunks = 1e3,
                   random_seed = NULL,
                   result_as_dt = FALSE){
    basicConfig(level = log.level)

    #check for NA values that make elasticnet crash
    if(anyNA(gene_expr))
        stop("NA values found in gene expression data. Can not proceed")
    if(anyNA(mir_expr))
        stop("NA values found in miRNA expression data. Can not proceed")

    if(is.null(mir_interactions)){
        logwarn("No information on miRNA gene interactions was provided,
                all miRNAs will be considered and runtime will likely explode.")
        mir_intersect <- colnames(mir_expr)

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

    #all pairwise combinations of selected genes
    if(is.null(gene.combinations)){
        loginfo("Computing all pairwise combinations of genes")
        #consider only genes that have miRNA interactions

        gene.combinations <-
           genes_pairwise_combinations(length(sel.genes))

        genes.as.indices <- TRUE
    }

    #make sure sel.genes order is used
    gene_expr <- gene_expr[,sel.genes]

    #for getting indices of things
    all_mirs <- colnames(mir_expr)
    all_genes <- colnames(gene_expr)

    loginfo("Beginning SPONGE run...")

    num_of_samples <- nrow(gene_expr)
    num_of_tasks <- min(nrow(gene.combinations), parallel.chunks)

    #convert to big.matrix for shared memory operations
    gene_expr_bm <- as.big.matrix(gene_expr)
    gene_expr_description <- describe(gene_expr_bm)

    mir_expr_bm <- as.big.matrix(mir_expr)
    mir_expr_description <- describe(mir_expr_bm)

    SPONGE_result <- foreach(gene_combis =
                    split_rows(gene.combinations,
                                            chunks = parallel.chunks),
                    i = iterators::icount(),
                    .combine=function(...) rbindlist(list(...)),
                    .multicombine=TRUE,
                    .packages = c("logging", "ppcor", "foreach", "iterators",
                                  "data.table", "bigmemory"),
                    .export = c("fn_get_shared_miRNAs"),
                    .options.RNG = random_seed) %dorng% {
        basicConfig(level = log.level)

        logdebug(paste("SPONGE: worker is processing chunk: ", i, sep=""))

        #attach bigmemory objects
        attached_gene_expr <- attach.big.matrix(gene_expr_description)
        attached_mir_expr <- attach.big.matrix(mir_expr_description)

        result <- foreach(gene_combi = iter(gene_combis, by="row"),
                          .combine=function(...) rbindlist(list(...))) %do% {

            if(genes.as.indices){
                geneA_idx <- gene_combi[1]
                geneB_idx <- gene_combi[2]
                geneA <- all_genes[geneA_idx]
                geneB <- all_genes[geneB_idx]
            }
            else {
                geneA <- as.character(gene_combi[1,1])
                geneB <- as.character(gene_combi[1,2])
                geneA_idx <- which(all_genes == geneA)
                geneB_idx <- which(all_genes == geneB)
            }

            logdebug(paste("Processing source gene", geneA,
                           "and target gene", geneB ))

            source_expr <- attached_gene_expr[,geneA_idx]
            target_expr <- attached_gene_expr[,geneB_idx]

            #check correlation
            dcor <- cor(source_expr, target_expr, use = "complete.obs")

            if(is.na(dcor)){
                logdebug(paste("Source gene", geneA, "or target gene", geneB,
                               "have a correlation of", dcor,
                               "there is no variance in one of the genes"))
                return(NULL)
            }
            if(!is.null(min.cor))
                if(dcor < min.cor){
                    logdebug(paste("Source gene", geneA, "and target gene", geneB,
                               "have a correlation of", dcor,
                               "which is below the threshold of", min.cor))
                    return(NULL)
            }

            #check if miRNA interaction information is provided, otherwise we
            #consider ALL miRNAs in each comparison
            if(!is.null(mir_interactions)){
                mir_intersect <- fn_get_shared_miRNAs(geneA, geneB,
                                                      mir_interactions)

                #check if shared miRNAs are in expression matrix
                if(length(setdiff(mir_intersect, all_mirs)) > 0){
                    logdebug(paste("Source gene", geneA, "and target gene", geneB,
                              "shared miRNAs not found in mir_expr are discarded"))
                    mir_intersect <- intersect(mir_intersect, all_mirs)
                }

                #check if there are actually any shared mirnas
                if(length(mir_intersect) == 0){
                    logdebug(paste("Source gene", geneA, "and target gene", geneB,
                            "do not share any miRNAs with significant regulation"))
                    return(NULL)
                }
            }

            if(each.miRNA){
                result <- foreach(mirna = mir_intersect,
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
        logdebug(paste("SPONGE finished chunk:", i, "of", num_of_tasks))

        return(as.data.table(result))
    }
    if(result_as_dt) return(SPONGE_result)
    else return(as.data.frame(SPONGE_result))
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
