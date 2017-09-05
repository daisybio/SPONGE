#generates a positive semi definite matrix
posdef <- function (n, ev = runif(n, 0, 10))
{
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
}

#computes the schur complement with respect to R22
schur <- function(R) {
    if(ncol(R) < 3) stop("Input matrix needs to have dimensions 3x3 or larger")
    if(nrow(R) != ncol(R)) stop("Input matrix needs to be square")
    ((R[1:2,1:2]) - (R[1:2,3:ncol(R)]) %*%
         (solve(R[3:nrow(R),3:ncol(R)])) %*% t(R[1:2,3:ncol(R)]))
}

#computes sensitivity correlation
get.q <- function(S) {
    if(ncol(S) < 3) stop("Input matrix needs to have dimensions 3x3 or larger")
    if(nrow(S) != ncol(R)) stop("Input matrix needs to be square")
    s11 <- S[1:2,1:2] #gene-gene covariance
    pva <- schur(S) #partial covariance
    q <- s11/outer(sqrt(diag(s11)),sqrt(diag(s11))) - pva /
        outer(sqrt(diag(pva)),sqrt(diag(pva))) #- sens - corr
    return(q)
}

#computes lambda from the dot product
checkLambda <- function(x,y){
    if(length(x) != length(y))
        stop("input should be two vectors of equal length")
    beta <- foreach(i = seq_along(x), .combine = sum) %do% {x[i] * y[i]}
    norm.x <- base::norm(x,type="2")
    norm.y <- base::norm(y,type="2")
    (beta / (norm.x * norm.y))
}

quadraticSolver <- function(a, b, c){
    bsquare4ac <- b^2 - 4 * a * c

    tryCatch({
    if(bsquare4ac < 0 || a == 0){
        return(NULL)
    }

    solution.1 <- (-b + sqrt(bsquare4ac)) / (2 * a)
    solution.2 <- (-b - sqrt(bsquare4ac)) / (2 * a)

    return(list(solution.1, solution.2))
    }, error = function(e){
        return(NULL)
    })
}

#' Sampling zero multiple miRNA sensitivity covariance matrices
#'
#' @param m number of miRNAs, i.e. number of columns of the matrix
#' @param number_of_solutions stop after this many instances have been samples
#' @param number_of_attempts give up after that many attempts
#' @param gene_gene_correlation optional, define the correlation of the first
#' two elements, i.e. the genes.
#' @param log.level the log level, typically set to INFO, set to DEBUG for
#' verbose logging
#' @importFrom gRbase cov2pcor
#' @import ppcor
#' @import expm
#' @import logging
#' @import foreach
#'
#' @return a list of covariance matrices with zero sensitivity correlation
#' @export
#'
#' @examples sample_zero_mscor_cov(m = 1,
#' number_of_solutions = 1,
#' gene_gene_correlation = 0.5)
sample_zero_mscor_cov <- function(m, number_of_solutions,
                                 number_of_attempts = 1e3,
                                 gene_gene_correlation = NULL,
                                 log.level = "OFF"){

    loginfo(
        paste0(
            "Sampling covariance matrices for cor = ",
            gene_gene_correlation
        )
    )
    solutions <- foreach(solution = seq_len(number_of_solutions),
                         .packages = c("MASS", "gRbase",
                                       "ppcor", "expm", "logging",
                                       "foreach"),
                         .export = c("quadraticSolver",
                                     "checkLambda",
                                     "get.q",
                                     "schur",
                                     "posdef")) %dopar% {
        total <- 0
        basicConfig(level = log.level)

        logdebug(
            paste("Looking for zero sensitivity covariance matrix for case m =",
                      m, "solution no.", solution))
        #a lot of solutions are not within our constraints, so we repeat
        while(total < number_of_attempts){
            total <- total + 1
            if(is.null(gene_gene_correlation)) K = runif(1,-1,1) # R11
            else K = gene_gene_correlation
            v1 <- runif(m, -1, 1) #random v1

            if(m == 1){
                R22 <- 1
                Z <- v1

                a <- K^2 + Z^2 - K^2 * Z^2
                b <- -2 * K * Z
                c <- Z^2 * K^2

                v2.solutions <- quadraticSolver(a, b, c)
                if(is.null(v2.solutions)){
                    logdebug("no solution for v2 found")
                    next
                }
                else{
                    v2 <- v2.solutions[[1]]
                }
            }
            else{
                lambda = runif(1,-1,1) #cos theta

                R22 <- cov2cor(posdef(m)) #miRNA correlation
                R22_invsqrt <- ginv(sqrtm(R22)) #R^-1/2

                u1 <- (R22_invsqrt %*% v1)[,1]
                u1.norm = base::norm(u1, type = "2")

                #solving quadratic equation to obtain u2 solutions
                a <- ((K^2 - lambda^2) * u1.norm^2 - K^2)
                b <- 2 * lambda * K * u1.norm
                c <- -K^2 * u1.norm^2

                u2.solutions <- quadraticSolver(a, b, c)
                if(is.null(u2.solutions)){
                    logdebug("no solution for u2")
                    next
                }

                u2.norm.1 <- u2.solutions[[1]]
                u2.norm.2 <- u2.solutions[[2]]

                r12.m <- (K - sqrt(u1.norm^2) * sqrt(u2.norm.1^2) * lambda) *
                    (1 - u1.norm^2)^-(1/2) * (1 - u2.norm.1^2)^-(1/2)

                if(is.nan(r12.m)){
                    r12.m <- (K - sqrt(u1.norm^2) * sqrt(u2.norm.2^2) *lambda) *
                    (1 - u1.norm^2)^-(1/2) * (1 - u2.norm.2^2)^-(1/2)
                }

                if(is.nan(r12.m)){
                    logdebug("||u2|| is not a valid solution")
                    next
                }

                if(r12.m != K){
                    logdebug("correlation is not equal to partial correlation for selected ||u2||")
                    next
                }

                constraints <- (R22_invsqrt %*% rep(1,m))[,1]

                if(anyNA(constraints)){
                    logdebug("R22^-(1/2) invalid. Can not compute constraints")
                    next
                }

                u1.m <- u1[m]
                u2.wo_m <- foreach(i = seq_len(m-1), .combine = c) %do%{
                    runif(1, -constraints[i], constraints[i])
                }

                if(anyNA(u2.wo_m)){
                    logdebug("R22^-(1/2) invalid. Can not compute constraints")
                    next
                }

                beta <- foreach(i = seq_len(m-1), .combine = sum) %do%
                    {u1[i] * u2.wo_m[i]}

                A <- u1.m^2 - lambda^2 * u1.norm^2
                B <- 2 * beta * u1.m
                C <- beta^2 - lambda^2  * u1.norm^2 * sum(u2.wo_m^2)

                u2.solutions <- quadraticSolver(A, B, C)

                if(is.null(u2.solutions)){
                    logdebug("no solution for the k-th element of u2")
                    next
                }
                #solve quadratic equation
                u2.k.1 <- u2.solutions[[1]]
                u2.k.2 <- u2.solutions[[2]]

                #check if lambda constraint is really fulfilled
                u2 <- c(u2.wo_m, u2.k.1)

                #if the sign is wrong we have to take the other solution
                if(checkLambda(u1, u2) != lambda)
                    u2 <- c(u2.wo_m, u2.k.2)

                #normalize and then multiply by ||u2||
                u2.scaled <- u2 / base::norm(u2, type = "2") * u2.norm.1

                within_constraints <- any(abs(u2.scaled) > constraints)

                if(!within_constraints){
                    logdebug("solution violates constraints")
                    next
                }
                else
                {
                    #compute v2
                    v2 = (solve(R22_invsqrt) %*% u2.scaled)[,1]
                }
            }

            #construct correlation matrix R
            R <- matrix(ncol = m+2, nrow = m+2)
            diag(R) <- 1
            R[1,2] <- K
            R[2,1] <- K
            R[3:(2+m),3:(2+m)] <- R22
            R[1,3:(2+m)] <- v1
            R[2,3:(2+m)] <- v2
            R[3:(2+m), 1] <- v1
            R[3:(2+m), 2] <- v2

            #generate covariance matrix S
            L = diag(exp(rnorm(2+m)))
            S = L %*% R %*% L

            #test for negative variance
            if(any(diag(schur(S)) < 0)){
                logdebug("negative variance in partial covariance matrix")
                next
            }

            #test sensitivity correlation is zero
            mscor <- get.q(S)[1,2]
            if(is.nan(mscor)){
                logdebug("sensitivity correlation is NaN")
                next
            }
            else if(abs(mscor) > sqrt(.Machine$double.eps)){
                logdebug("sensitivity correlation is not zero")
                next
            }
            else{
                loginfo(paste("viable solution found for m =", m,
                              "and k =", K,
                               "solution no.", solution))
                attr(S, "iterations") <- total
                attr(S, "k") <- K
                attr(S, "m") <- m
                return(S)
            }
        }
    }
    solutions <- Filter(Negate(is.null), solutions)
    if(length(solutions) == 0){
        logerror("No solutions found")
        return(NULL)
    }
    else{
        total <- sum(unlist(lapply(solutions, function(x){
            attr(x, "iterations")})))
        loginfo(paste("case k = ", gene_gene_correlation, "m = ", m, " - Found",
                      length(solutions), "solutions in a total of",
                      total, "iterations."))
    }

    return(solutions)
}


#' Sample mscor coefficients from pre-computed covariance matrices
#'
#' @param cov_matrices a list of pre-computed covariance matrices
#' @param number_of_samples  the number of samples available in the expression
#' data
#' @param number_of_datasets the number of mscor coefficients to be sampled
#' from each covariance matrix
#' @seealso sample_zero_mscor_cov
#' @return a vector of mscor coefficients
#' @export
#' @import MASS
#' @import foreach
#' @import gRbase
#' @import ppcor
#'
#' @examples #we select from the pre-computed covariance matrices in SPONGE
#' #10 for m = 5 miRNAs and gene-gene correlation 0.6
#' cov_matrices_selected <- precomputed_cov_matrices[["5"]][["0.6"]]
#' sample_zero_mscor_data(cov_matrices = cov_matrices_selected,
#' number_of_samples = 200, number_of_datasets = 10)
sample_zero_mscor_data <- function(cov_matrices,
                                  number_of_samples = 100,
                                  number_of_datasets = 100){
    foreach(cov.matrix = cov_matrices) %do% {
        #check that sensitivity correlation is zero
        if(abs(cov2pcor(cov.matrix)[1,2] - cov2cor(cov.matrix)[1,2]) >
           sqrt(.Machine$double.eps))
            stop("sensitivity correlation of a given covariance matrix is not zero.")

        #sample data under this covariance matrix
        foreach(i = seq_len(number_of_datasets), .combine = c) %do%{
            sample.data <- mvrnorm(n = number_of_samples,
                                   rep(0, ncol(cov.matrix)),
                                   cov.matrix,
                                   empirical = FALSE)
            cor(sample.data)[1,2] - pcor(sample.data)$estimate[1,2]
        }
    }
}

