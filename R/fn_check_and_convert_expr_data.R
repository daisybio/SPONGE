#' Checks if expression data is in matrix or ExpressionSet format and
#' converts the latter to a standard matrix. Alternatively, a big.matrix
#' descriptor object can be supplied to make use of shared memory between
#' parallelized workers through the bigmemory package.
#'
#' @param expr_data expr_data as matrix or ExpressionSet
#' @importFrom Biobase exprs
#' @return expr_data as matrix
#'
#' @examples \dontrun{check_and_convert_expression_data(gene_expr)}
check_and_convert_expression_data <- function(expr_data){

    if(is(expr_data, "big.matrix.descriptor")){
        expr_data <- attach.big.matrix(expr_data)
        if(length(mwhich(expr_data, 1:ncol(expr_data), NA, 'eq', 'OR')) > 0){
            stop("NA values found in expression data. Can not proceed")
        }
        else return(expr_data)
    }

    if(!is(expr_data, "matrix") && !is(expr_data, "ExpressionSet")){
        stop("expression matrix has to be of class matrix, ExpressionSet or big.matrix.descriptor")
    }

    if(is(expr_data, "ExpressionSet")){
        expr_data <- t(exprs(expr_data))
    }

    #check for NA values that make elasticnet crash
    if(anyNA(expr_data)) stop("NA values found in expression data. Can not proceed")

    return(expr_data)
}