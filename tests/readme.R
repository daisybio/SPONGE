#we use testthat, see subfolder for unit tests.
if(!(requireNamespace("testthat", quietly = TRUE)))
    stop("install package testthat to run tests")

if(!(requireNamespace("doParallel", quietly = TRUE)))
    stop("install package doParallel to run tests")