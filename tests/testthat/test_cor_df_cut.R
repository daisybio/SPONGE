library(SPONGE)

context("Test SPONGE assigning ceRNA interactions to partitions correctly")

test_that("test if cor_cut and df_cut are selecte appropriately",{
    set.seed(12345)
    partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = seq(0.2, 0.90, 0.1), #gene-gene correlation
        m_max = 8, #number of miRNAs
        sponge_result = ceRNA_interactions)

    #expect that correlations are in a +-0.05 window around their cor_cut
    expect_false(any(abs(as.numeric(as.character(partitions$cor_cut)) -
                                    partitions$cor) > 0.1))
    partitions_df_up_to_8 <- partitions[df < 9]
    partitions_df_more_than_8 <- partitions[df > 8]

    #expect that up to a df of 8, df_cut corresponds to df
    expect_equal(as.numeric(as.character(partitions_df_up_to_8$df_cut)),
                        partitions_df_up_to_8$df)

    #for larger df use the model for df = 8
    expect_false(any(as.numeric(as.character(partitions_df_more_than_8$df_cut)) >
                        partitions_df_more_than_8$df))
})

test_that("test if cor_cut and df_cut are selecte appropriately with large m",{
    set.seed(12345)
    partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = seq(0.2, 0.90, 0.1), #gene-gene correlation
        m_max = 15, #number of miRNAs
        sponge_result = ceRNA_interactions)

    #expect that correlations are in a +-0.05 window around their cor_cut
    expect_false(any(abs(as.numeric(as.character(partitions$cor_cut)) -
                             partitions$cor) > 0.1))
    partitions_df_up_to_8 <- partitions[df < 9]
    partitions_df_more_than_8 <- partitions[df > 8]

    #expect that up to a df of 8, df_cut corresponds to df
    expect_equal(as.numeric(as.character(partitions_df_up_to_8$df_cut)),
                 partitions_df_up_to_8$df)

    #for larger df use the model for df = 8
    expect_false(any(as.numeric(as.character(partitions_df_more_than_8$df_cut)) >
                         partitions_df_more_than_8$df))
})


test_that("test if cor_cut and df_cut are selecte appropriately with custom k",{
    set.seed(12345)
    partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = c(0.2, 0.5, 0.7), #gene-gene correlation
        m_max = 8, #number of miRNAs
        sponge_result = ceRNA_interactions)

    #expect that correlations are in a +-0.05 window around their cor_cut
    expect_false(any(abs(as.numeric(as.character(partitions$cor_cut)) -
                             partitions$cor) > 0.15))
    partitions_df_up_to_8 <- partitions[df < 9]
    partitions_df_more_than_8 <- partitions[df > 8]

    #expect that up to a df of 8, df_cut corresponds to df
    expect_equal(as.numeric(as.character(partitions_df_up_to_8$df_cut)),
                 partitions_df_up_to_8$df)

    #for larger df use the model for df = 8
    expect_false(any(as.numeric(as.character(partitions_df_more_than_8$df_cut)) >
                         partitions_df_more_than_8$df))
})

test_that("test if computing cor_cut and df_cut fails with k = 0",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = seq(0, 0.90, 0.1), #gene-gene correlation
        m_max = 15, #number of miRNAs
        sponge_result = ceRNA_interactions))
})

test_that("test if computing cor_cut and df_cut fails with k = 1",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = seq(0.1, 1, 0.1), #gene-gene correlation
        m_max = 15, #number of miRNAs
        sponge_result = ceRNA_interactions))
})


test_that("test if computing cor_cut and df_cut fails with k = NA",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = NA, #gene-gene correlation
        m_max = 15, #number of miRNAs
        sponge_result = ceRNA_interactions))
})

test_that("test if computing cor_cut and df_cut fails with m = NA",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = c(0.3, 0.5, 0.7), #gene-gene correlation
        m_max = NA, #number of miRNAs
        sponge_result = ceRNA_interactions))
})


test_that("test if computing cor_cut and df_cut fails with k = NULL",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = NULL, #gene-gene correlation
        m_max = 15, #number of miRNAs
        sponge_result = ceRNA_interactions))
})

test_that("test if computing cor_cut and df_cut fails with m = NULL",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = c(0.3, 0.5, 0.7), #gene-gene correlation
        m_max = NULL, #number of miRNAs
        sponge_result = ceRNA_interactions))
})

test_that("test if computing cor_cut and df_cut fails with missing columns",{
    set.seed(12345)
    expect_error(partitions <- determine_cutoffs_for_null_model_partitioning(
        ks = c(0.3, 0.5, 0.7), #gene-gene correlation
        m_max = NULL, #number of miRNAs
        sponge_result = data.frame(mscor = rnorm(5))))
})