Sys.unsetenv("R_TESTS")

library(SPONGE)
library(doParallel)
cl <- makeCluster(2)

context("TEST functions for generating cov matrices under the NULL")

test_that("Test that posdef returns positive semi-definite matrices",{
    for(i in 3:10){
        m <- posdef(i)
        expect_equal(ncol(m), i)
        expect_equal(nrow(m), i)
        expect_true(all(eigen(m)$values > 0))

    }
})

test_that("Schur complement works for a simple case", {
    m <- matrix(c(1,6,7,6,1,8,7,8,1),3,3)
    expect_equal(as.vector(schur(m)), c(-48, -50, -50, -63))
})

test_that("Computing sensitivity correlation works as expected", {

    for(i in seq_len(8)){
        mscor <- get.q(precomputed_cov_matrices[[i]][[sample(seq_len(8), 1)]][[1]])[1,1]
        expect_equal(mscor, 0)
    }
})

test_that("Check lambda function works as expected",{

    expect_equal(checkLambda(x = c(1,2,3), y = c(1,2,3)), 1)
    expect_equal(checkLambda(x = c(1,1,1), y = c(2,2,2)), 1)
    expect_equal(checkLambda(x = c(1,0), y = c(0,1)), 0)
})

test_that("Solving quadrativ equations works as expected", {

    expect_equal(unlist(quadraticSolver(a = 1, b = 2, c = -3)), c(1, -3))
})

test_that("Sampling cov. matrices works for m = 1", {

    simple_case <- sample_zero_mscor_cov(m = 1,number_of_solutions = 1,
                          gene_gene_correlation = 0.5,
                          random_seed = 12345)[[1]]
    expect_equal(as.vector(simple_case)[1], 1.409278, tolerance = 1e-6)
    expect_equal(ncol(simple_case), 3)

})

test_that("Sampling cov. matrices works for m = 3", {

    complex_case <- sample_zero_mscor_cov(m = 3,number_of_solutions = 1,
                                         gene_gene_correlation = 0.5,
                                         random_seed = 1274)[[1]]
    expect_equal(complex_case[1,1], 2.3165434, tolerance = 1e-6)
    expect_equal(ncol(complex_case), 5)
})

test_that("Sampling cov. matrices works for m = 3 in parallel", {

    registerDoParallel(cl)

    complex_case <- sample_zero_mscor_cov(m = 3, number_of_solutions = 1,
                                              gene_gene_correlation = 0.5,
                                              random_seed = 1274)[[1]]
    expect_equal(complex_case[1,1], 2.3165434, tolerance = 1e-6)
    expect_equal(ncol(complex_case), 5)

    registerDoSEQ()
})

test_that("Sampling data from covariance matrices works", {
    set.seed(12345)
    result <- unlist(
        sample_zero_mscor_data(
            cov_matrices = precomputed_cov_matrices[[1]][[1]][1:2],
            number_of_samples = 50,
            number_of_datasets = 2))
    expect_equal(result, c(0.05027585, 0.06882558, -0.05860923, -0.00840041))
})

stopCluster(cl)