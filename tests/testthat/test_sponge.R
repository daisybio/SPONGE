library(SPONGE)

context("TEST sponge prediction method")

test_that("sponge with mir interactions works",{

    result <- sponge(gene_expr = gene_expr[,1:3],
       mir_expr = mir_expr,
       mir_interactions = mir_interactions,
       random_seed = 1234)

    expect_equal(nrow(result), 3)
    expect_equal(result$geneA, c("UST", "UST", "FBXO32"))
    expect_equal(result$mscor, c(0.2275383, 0.0972899, 0.1855703),
                 tolerance = 1e-7)
})

test_that("sponge correlation filter works",{

    result <- sponge(gene_expr = gene_expr[,1:3],
                     mir_expr = mir_expr,
                     mir_interactions = mir_interactions,
                     min.cor = 0.2,
                     random_seed = 1234)

    expect_equal(nrow(result), 1)
    expect_equal(result$geneA, c("UST"))
    expect_equal(result$mscor, c(0.2275383),
                 tolerance = 1e-7)
})

test_that("sponge without correlation filter works",{

    result <- sponge(gene_expr = gene_expr[,1:3],
                     mir_expr = mir_expr,
                     mir_interactions = mir_interactions,
                     min.cor = NULL,
                     random_seed = 1234)

    expect_equal(nrow(result), 3)
})

test_that("sponge without mir interactions works",{

    result <- sponge(gene_expr = gene_expr[,1:3],
                     mir_expr = mir_expr,
                     mir_interactions = NULL,
                     random_seed = 1234)

    expect_equal(nrow(result), 3)
    expect_equal(result$geneA, c("UST", "UST", "FBXO32"))
    expect_equal(result$df, rep(50,3))
    expect_true(all(result$mscor - c(0.2275383, 0.0972899, 0.1855703) < 0))
})

test_that("sponge for indiviudal mirnas works",{

    result <- sponge(gene_expr = gene_expr[,1:3],
                     mir_expr = mir_expr,
                     mir_interactions = mir_interactions,
                     each.miRNA = TRUE,
                     random_seed = 1234)

    expect_equal(nrow(result), 12)
    expect_equal(result$df, rep(1,12))
    expect_equal(ncol(result), 7)
})

test_that("manually selecting genes works",{
    expect_equal(
        sponge(gene_expr = gene_expr,
               mir_expr = mir_expr,
               mir_interactions = mir_interactions,
               selected.genes = c("TMEM132A", "ESR1", "DNMT3A"),
               random_seed = 1234)$mscor,
        0.1138124,
        tolerance = 1e-7
    )
})

test_that("manually defining gene combinations works",{
    expect_equal(
        sponge(gene_expr = gene_expr,
               mir_expr = mir_expr,
               mir_interactions = mir_interactions,
               gene.combinations = data.frame(geneA = "DNMT3A",
                                              geneB = "TMEM132A"),
               random_seed = 1234)$mscor,
        0.1138124,
        tolerance = 1e-7
    )
    expect_equal(
        sponge(gene_expr = gene_expr,
               mir_expr = mir_expr,
               mir_interactions = mir_interactions,
               gene.combinations = data.frame(geneA = "TMEM132A",
                                              geneB = "DNMT3A"),
               random_seed = 1234)$mscor,
        0.1138124,
        tolerance = 1e-7
    )
})

test_that("combinations with empty results work",{
    expect_equal(
        nrow(sponge(gene_expr = gene_expr,
               mir_expr = mir_expr,
               mir_interactions = mir_interactions,
               gene.combinations = data.frame(geneA = "ESR1",
                                              geneB = "DNMT3A"),
               random_seed = 1234)), 0)
})

test_that("split rows works", {
  iterator <- split_rows(ceRNA_interactions, chunks = 4)
  firstElt <- nextElem(iterator)
  expect_equal(nrow(firstElt), 48)
  for(i in seq_len(3)) nextElem(iterator)
  expect_error(nextElem(iterator))
})

test_that("computing partial correlations works",{
    dcor <- cor(gene_expr[,1], gene_expr[,2])
    pcor <- pcor.test(gene_expr[,1], gene_expr[,2], mir_expr[,1])$estimate
    expect_equal(compute_pcor(source_expr = gene_expr[,1],
                 target_expr = gene_expr[,2],
                 m_expr = mir_expr[,1],
                 geneA = "A",
                 geneB = "B",
                 dcor = dcor)$mscor,
                 dcor - pcor)
})

test_that("mapping mimats to mir names works",{
    stop("implement me")
})

test_that("getting shared miRNAs of two genes works",{
    stop("implement me")
})