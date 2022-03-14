#' Compute the residual sum of squares error for an elastic net model
#' @param model The elastic net model
#' @param x The miRNA expression
#' @param y The gene expression
#'
#' @return the RSS
fn_get_rss <- function(model, x, y){
    lambda.min <- model$lambda.min
    #get coefficients (beta) of best model
    model.pred <- predict(model, x, s=lambda.min)
    model.rss <- sum((model.pred - y)^2)
    model.rss <- model.rss / ncol(x)
    return(model.rss)
}

#' Extract the model coefficients from an elastic net model
#' @param model An elastic net model
#' @import logging
#'
#' @return A data frame with miRNAs and coefficients
fn_get_model_coef <- function(model){
    lambda.min <- model$lambda.min
    model.coef <- coef(model, s=lambda.min)

    #extract names of those miRNAs where the coefficient is != 0
    #and remove intercept [-1]
    logdebug(paste("Extracting miRNAs with non zero coefficient for gene", gene))

    coefficients <- as.vector(model.coef)[-1]
    mimats <- rownames(model.coef)[-1]
    non.zero.model.coef <- which(coefficients != 0)

    #check if any of the coefficients were non-zero
    if(length(non.zero.model.coef)==0){
        logwarn("no non-zero coefficients for this model")
        return(NULL)
    }
    else{
        coefficients <- coefficients[non.zero.model.coef]
        mimats <- mimats[non.zero.model.coef]
        data.frame(mirna=mimats, coefficient=coefficients)
    }
}

#' Computes an elastic net model
#' @import foreach
#' @importFrom glmnet cv.glmnet
#' @param x miRNA expression matrix
#' @param y gene expression vector
#' @param alpha.step Step size for alpha, the tuning parameter for elastic net.
#'
#' @return The best model, i.e. the one for which the selected alpha yielded the
#' smallest residual sum of squares error
fn_elasticnet <- function(x, y, alpha.step = 0.1){
    models <- foreach(alpha = seq(0, 1, alpha.step)) %do%{
        tryCatch({
            glmnet::cv.glmnet(x, y, alpha = alpha)
        }, warning = function(w){
            logwarn(w)
            return(glmnet::cv.glmnet(x, y, alpha = alpha))
        }, error = function(e){
            logerror(e)
            return(NA)
        })
    }

    models.cvm <- unlist(lapply(models, function(model){
        min(model$cvm)
    }))

    #return model with smallest residual sum of squares
    return(models[[which.min(models.cvm)]])
}

