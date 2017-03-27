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
#' @importFrom matrixStats colVars
#' @import targetscan.Hs.eg.db
#' @import dplyr
#' @import foreach
#' @import glmnet
#' @import bigmemory
#'
#' @param gene_expr A matrix of gene expression
#' @param mir_expr A matrix of miRNA expression
#' @param mir_predicted_targets A data frame with miRNA in cols and genes in rows.
#' A 0 indicates the miRNA is not predicted to target the gene, >0 otherwise.
#' If this parameter is NULL all miRNA-gene interactions are tested
#' @param Some miRNA target databases do not use mature miRNAs but miRNA family
#' ids since members of the same family have an identical seed sequence. To map
#' miRNA expression to such a database we need a mapping table. It can be
#' created by downloading the miR family info from the mirbase db and converting
#' it with the function build_mir_family_dictionary
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
#' @param coefficient.threshold If F.test is FALSE, the
#' regression model considers only negative coefficients < threshold.
#' This is a sensible strategy since only miRNA with a negative coefficient
#' have an inhibitory role on gene expression.
#' @return A list of genes, where for each gene, the regulating miRNA are
#' included as a data frame. For F.test = TRUE this is a data frame with fstat
#' and p-value for each miRNA. Else it is a data frame with the model
#' coefficients.
#' @export
#'
#' @seealso build_mir_family_dictionary
#' @seealso sponge
#' @examples
#' #library(doParallel)
#' #cl <- makePSOCKcluster(2)
#' #registerDoParallel(cl)
#' genes_miRNA_candidates <- gene_miRNA_interaction_filter(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' mir_predicted_targets = NULL)
#' #stopCluster(cl)
#'
#' #If we also perform an F-test, only few of the above miRNAs remain
#' genes_miRNA_candidates <- gene_miRNA_interaction_filter(
#' gene_expr = gene_expr,
#' mir_expr = mir_expr,
#' F.test = TRUE,
#' F.test.p.adj.threshold = 0.05,
#' mir_predicted_targets = NULL)
#'
gene_miRNA_interaction_filter <- function(gene_expr, mir_expr,
                                         mir_predicted_targets = list(mircode, targetscan),
                                         elastic.net = TRUE,
                                         log.every.n = 10,
                                         log.level='INFO',
                                         var.threshold = NULL,
                                         F.test = FALSE,
                                         F.test.p.adj.threshold = 0.05,
                                         coefficient.threshold = 0,
                                         randomize.predicted.targets = FALSE){

    #remove columns with little variance
    if(!is.null(var.threshold)){
        gene_expr <- gene_expr[,which(colVars(gene_expr) > var.threshold)]
        mir_expr <- mir_expr[,which(colVars(mir_expr) > var.threshold)]
    }
    num_of_genes <- ncol(gene_expr)

    #convert to big.matrix for shared memory operations
    gene_expr <- as.big.matrix(gene_expr)
    gene_expr_description <- describe(gene_expr)

    mir_expr <- as.big.matrix(mir_expr)
    mir_expr_description <- describe(mir_expr)

    #merge mirna target annotation
    loginfo("merging miRNA target database annotations")
    if(is.list(mir_predicted_targets)){
        list_of_mir_predicted_targets <- mir_predicted_targets

        all_mirs <- foreach(mir_db = mir_predicted_targets,
                            .combine = union) %do% {
                                colnames(mir_db)
                            }

        all_genes <- foreach(mir_db = mir_predicted_targets,
                            .combine = union) %do% {
                                rownames(mir_db)
                            }

        mir_predicted_targets <- big.matrix(nrow = length(all_genes),
                             ncol = length(all_mirs),
                             dimnames = list(all_genes, all_mirs))

        mir_predicted_targets[,] <- 0

        #fill big matrix
        for(mir_db in list_of_mir_predicted_targets){
            mir_predicted_targets[rownames(mir_db), colnames(mir_db)] <-
                mir_predicted_targets[rownames(mir_db), colnames(mir_db)] + mir_db
        }
    }
    else{
        mir_predicted_targets <- as.big.matrix(mir_predicted_targets)
    }

    #randomize target genes
    if(randomize.predicted.targets){
        rownames(mir_predicted_targets) <- sample(rownames(mir_predicted_targets))
    }

    mir_predicted_targets_description <- describe(mir_predicted_targets)

    #loop over all genes and compute regression models to identify important miRNAs
    foreach(col.num = 1:ncol(gene_expr),
            .final = function(x) setNames(x, colnames(gene_expr)),
            .packages = c("logging", "glmnet", "dplyr", "bigmemory"),
            .export = c("fn_get_model_coef", "fn_elasticnet", "fn_gene_miRNA_F_test", "fn_get_rss"),
            .inorder = TRUE) %dopar% {

        #attach shared data
        attached_gene_expr <- attach.big.matrix(gene_expr_description)
        attached_mir_expr <- attach.big.matrix(mir_expr_description)
        attached_mir_predicted_targets <- attach.big.matrix(mir_predicted_targets_description)

        #get gene name
        gene <- colnames(attached_gene_expr)[col.num]

        #setup logging
        basicConfig(level = log.level)

        curr_percentage <- round((col.num / num_of_genes) * 100, 2)
        if(col.num %% log.every.n == 0) loginfo(paste("Computing gene / miRNA regression models: ", curr_percentage, "% completed.", sep=""))

        logdebug(paste("Processing gene", gene))

        #expression values of this gene
        g_expr <- attached_gene_expr[,gene]

        #if there is zero variance this exercise is pointless
        if(var(g_expr) == 0){
            logdebug(paste("Zero variance found in gene", gene, ". Returning null for this gene"))
            return(NULL)
        }

        #check if we have a target database
        if(!is.null(attached_mir_predicted_targets)){
            #miRNAs this genes has binding sites for
            if(!(gene %in% rownames(attached_mir_predicted_targets))){
                logdebug(paste("No information about miRNA target interaction found for ", gene, ". Returning null for this gene", sep=""))
                return(NULL)
            }
            else{
                logdebug(paste("Extracting miRNAs targeting gene", gene))
            }

            mimats <- colnames(attached_mir_predicted_targets)[which(attached_mir_predicted_targets[gene,] > 0)]

            if(length(mimats) == 0){
                logdebug("Gene not found in miRNA target database. Returning null for this gene.")
                return(NULL)
            }

            #expression of miRNAs selected above
            mimats_matched <- intersect(mimats, colnames(attached_mir_expr))

            if(length(mimats_matched) == 0){
                logdebug("None of the target mirnas are found in expression data. Returning null for this gene.")
                return(NULL)
            }
            else if(!elastic.net){
                return(data.frame(mirna = mimats_matched))
            }
            logdebug(paste(length(mimats_matched), "of", length(mimats), "mirnas found in expression data for gene", gene))

            m_expr <- attached_mir_expr[,mimats_matched]
        }
        else{
            logwarn("No miRNA target database provided. Using all miRNAs")
            mimats_matched <- colnames(attached_mir_expr)
            m_expr <- attached_mir_expr
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
                    lm_result <- data.frame(mirna=mimats_matched[1], fstats=fstat, pval=pval, p.adj = pval)
                }else{
                    lm_result <- data.frame(mirna = mimats_matched[1], coefficient = coef(lm(g_expr ~ m_expr))[-1])
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
            smaller.than.threshold <- which(result$coefficient < coefficient.threshold)
            if(length(smaller.than.threshold) == 0) return(NULL)
            else return(result[smaller.than.threshold,])
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
