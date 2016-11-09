#' Perform F test for gene-miRNA elastic net model
#'
#' @param gene_expr A gene expression matrix with samples in rows and
#' genes in columns
#' @param mir_expr A miRNA expression matrix with samples in rows and
#' genes in columns. Sample number and order has to agree with
#' above gene expression matrix
#' @param model A nested elastic net model to be tested
#'
#' @return return data frame with miRNA, fstat and adjusted p.value (BH).
#'
fn_gene_miRNA_F_test <- function(gene_expr, mir_expr, model,
                                 p.adj.threshold = NULL){

    model.rss <- fn_get_rss(model, mir_expr, gene_expr)

    non.zero.model.coef.mirnas <- fn_get_model_coef(model)[,1]

    #remove each of the relevant features to assess its signifiance in an F test
    model.stats <- foreach(mirna = non.zero.model.coef.mirnas, .combine=rbind) %do%{
        #set this feature to zero to ignore it in the model
        test_mir_expr <- mir_expr
        test_mir_expr[,mirna] <- 0

        # use RSS to compute the F-statistics for nested model.
        test.rss <- fn_get_rss(model, test_mir_expr, gene_expr)

        # F test
        df2 <- nrow(mir_expr) - length(non.zero.model.coef.mirnas) - 1 #samples - features

        #avoid negative dof
        if(df2 < 0) df2 <- 0
        fstats <- (test.rss - model.rss)/(model.rss/df2)

        data.frame(mirna=mirna,
                   fstats=fstats,
                   pval=pf(fstats, 1, df2, lower.tail=FALSE))
    }

    if(!is.null(model.stats)){
        model.stats$p.adj <- p.adjust(model.stats$pval, "BH")
        if(!is.null(p.adj.threshold)){
            model.stats <- dplyr::filter(model.stats, p.adj < p.adj.threshold)
        }
    }

    return(model.stats)
}