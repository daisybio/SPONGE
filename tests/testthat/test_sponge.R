library(SPONGE)

context("Test sponge prediction method")

test_that("manually defining gene combinations works",{
          expect_equal(
              sponge(gene_expr = gene_expr,
                     mir_expr = mir_expr,
                     mir_interactions = mir_interactions,
                     gene.combinations = data.frame(geneA = "CDH1",
                                                    geneB = "SERINC1"))$pcor,
              0.175298,
              tolerance = 1e-7
          )
})