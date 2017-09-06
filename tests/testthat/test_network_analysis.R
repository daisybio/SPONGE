library(SPONGE)

context("TEST network analysis functions")

test_that("node centralities can be computed",{
    ncs <- sponge_node_centralities(sponge_result = ceRNA_interactions)

    expect_equal(nrow(ncs), 25)
    expect_equal(ncol(ncs), 5)
    expect_true(is.numeric(ncs$degree))
    expect_true(is.numeric(ncs$eigenvector))
    expect_equal(max(ncs$eigenvector), 1)
    expect_true(is.numeric(ncs$betweenness))
    expect_true(is.numeric(ncs$page_rank))
})

test_that("edge centralities can be computed",{
    ecs <- sponge_edge_centralities(sponge_result = ceRNA_interactions)

    expect_equal(nrow(ecs), 191)
    expect_equal(ncol(ecs), 3)
    expect_true(is.numeric(ecs$edge_betweenness))
})