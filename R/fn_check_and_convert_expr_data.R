#' Checks if expression data is in matrix or ExpressionSet format and
#' converts the latter to a standard matrix.
#'
#' @param expr_data expr_data as matrix or ExpressionSet
#' @importFrom Biobase exprs
#' @return expr_data as matrix
#'
#' @examples check_and_convert_expression_data(gene_expr)
check_and_convert_expression_data <- function(expr_data){
    if(!class(expr_data) %in% c("matrix", "ExpressionSet"))
        stop("expression matrix has to be of class matrix or ExpressionSet")

    if(class(expr_data) == "ExpressionSet"){
        return(t(exprs(expr_data)))
    }
    else{
        return(expr_data)
    }
}