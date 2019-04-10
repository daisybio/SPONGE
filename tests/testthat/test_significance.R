library(SPONGE)
library(iterators)

context("TEST significance functions")

test_that("drawing sample data from a set of covariance matrices works",{
    set.seed(12345)
    selected_cov_matrices <- precomputed_cov_matrices[[1]]
    expect_equal(
        compute_null_model(cov_matrices = selected_cov_matrices[[1]],
                               number_of_datasets = 5,
                               number_of_samples = 100)$mscor,
        c(-0.1629024561, -0.1410354380, -0.0589320903, 0.0008393662, 0.0812918164),
        tolerance = 1e-7
    )
})

test_that("test computing a null model", {
    null_model <- sponge_build_null_model(
        number_of_datasets = 5,
        m_max = 3,
        number_of_samples = 20,
        cov_matrices = precomputed_cov_matrices[1:3])

    expect_length(null_model, 3)
    expect_true(all(unlist(lapply(null_model, length)) == 8))
})

test_that("computing a null model fails when sample number too small", {
    expect_error(sponge_build_null_model(number_of_datasets = 5,
                                          number_of_samples = 10,
                                          cov_matrices = precomputed_cov_matrices))
})

set.seed(12345)

null_model <- sponge_build_null_model(
    number_of_datasets = 5,
    m_max = 3,
    number_of_samples = 20,
    cov_matrices = precomputed_cov_matrices[1:3])

test_that("isplitDT2 returns a working iterator",{

    ks <- seq(0.2, 0.90, 0.1)

    partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = ks,
        m_max = 3,
        sponge_result = ceRNA_interactions)

    iterator <- isplitDT2(x = partitions, ks = ks, ms = seq(1,3,1),
                          null_model = null_model)
    first_elt <- nextElem(iterator)

    expect_equal(first_elt$key, c(Var1 = 0.2, Var2 = 1.0))
    expect_true(nrow(first_elt$sim.data) == 5)
    expect_true(nrow(first_elt$value) == 10)
    expect_true(all(first_elt$value$cor_cut == "0.2"))
    expect_true(all(first_elt$value$df_cut == "1"))
    for(i in seq_len(((length(ks) * 3) - 1))) nextElem(iterator)
    expect_error(nextElem(iterator))
})

test_that("isplitDT2 throws error if simulated data is missing",{

    ks <- seq(0.2, 0.90, 0.1)

    partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = ks,
        m_max = 3,
        sponge_result = ceRNA_interactions)

    iterator <- isplitDT2(x = partitions, ks = ks, ms = seq(1,8,1),
                          null_model = null_model)
    for(i in seq_len(8*3)) nextElem(iterator)
    expect_error(nextElem(iterator))
})

test_that("computing p-values is working", {
    p_vals <- sponge_compute_p_values(sponge_result = ceRNA_interactions,
                                      null_model = null_model)
    expect_true(all(p_vals$p.val <= 1))
    expect_true(all(p_vals$p.val > 0))
    expect_true(all(p_vals$p.adj >= p_vals$p.val))
})
