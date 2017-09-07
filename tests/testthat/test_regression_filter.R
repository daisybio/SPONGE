library(SPONGE)
library(doParallel)
cl <- makeCluster(2)

context("TEST gene miRNA regression filter")

test_that("Regression filter works with one interaction db",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = targetscan_symbol,
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(genes_miRNA_candidates$ASAP2$coefficient, -0.2553748,
                 tolerance = 1e-7)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 7)

    #test if interactions are really in mircode
    expect_true(all(as.character(genes_miRNA_candidates$ASAP2$mirna) %in%
                        names(which(mircode_symbol["ASAP2",] > 0))))
    expect_true(all(as.character(genes_miRNA_candidates$DNMT3A$mirna) %in%
                        names(which(mircode_symbol["DNMT3A",] > 0))))
})

test_that("Regression filter accepts ExpressionSet as input",{
    gene_expr_set <- Biobase::ExpressionSet(assayData =
                                    t(gene_expr[,c("ASAP2", "DNMT3A")]))
    mir_expr_set <- Biobase::ExpressionSet(assayData = t(mir_expr))

    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr_set,
        mir_expr = mir_expr_set,
        mir_predicted_targets = targetscan_symbol,
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(genes_miRNA_candidates$ASAP2$coefficient, -0.2553748,
                 tolerance = 1e-7)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 7)

    #test if interactions are really in mircode
    expect_true(all(as.character(genes_miRNA_candidates$ASAP2$mirna) %in%
                        names(which(mircode_symbol["ASAP2",] > 0))))
    expect_true(all(as.character(genes_miRNA_candidates$DNMT3A$mirna) %in%
                        names(which(mircode_symbol["DNMT3A",] > 0))))
})

test_that("Regression filter works with stricter threshold",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = targetscan_symbol,
        coefficient.threshold = -0.2,
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(genes_miRNA_candidates$ASAP2$coefficient, -0.2553748,
                 tolerance = 1e-7)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 2)
})

test_that("Regression filter no results with > sign for threshold",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = targetscan_symbol,
        coefficient.threshold = 0,
        coefficient.direction = ">",
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_true(all(genes_miRNA_candidates$DNMT3A$coefficient >0))
})

test_that("Regression filter no results with > sign and cutoff for threshold",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = targetscan_symbol,
        coefficient.threshold = 0.2,
        coefficient.direction = ">",
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 1)
})

test_that("Regression filter works without mir interaction db", {
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = NULL,
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 18)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 7)
})

test_that("Regression filter works with several mir interaction dbs",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = list(mircode_symbol, targetscan_symbol),
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 16)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 6)
})

test_that("Regression filter works without threshold",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = list(mircode_symbol, targetscan_symbol),
        coefficient.threshold = 0,
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 25)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 13)
})

test_that("Regression filter works with F-test",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        F.test = TRUE,
        F.test.p.adj.threshold = 1,
        mir_predicted_targets = list(mircode_symbol, targetscan_symbol),
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 36)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 12)
})

test_that("Regression filter works with F-test and p.adj cutoff",{
    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        F.test = TRUE,
        F.test.p.adj.threshold = 1e-50,
        mir_predicted_targets = list(mircode_symbol, targetscan_symbol),
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 13)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 7)
})


test_that("Regression filter works in no targets mode",{

    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = list(mircode_symbol),
        random_seed = 1234,
        select.non.targets = TRUE)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 3)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 5)

    #check if selected miRNAs are really not reported as targets in mircode
    expect_true(all(as.character(genes_miRNA_candidates$ASAP2$mirna) %in%
                        names(which(mircode_symbol["ASAP2",] == 0))))
    expect_true(all(as.character(genes_miRNA_candidates$DNMT3A$mirna) %in%
                        names(which(mircode_symbol["DNMT3A",] == 0))))
})

test_that("Regression filter with variance threshold",{

    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        var.threshold = 0.5,
        mir_predicted_targets = list(mircode_symbol),
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 12)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 6)
})

test_that("Regression filter with too strict variance threshold",{

    expect_error(genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        var.threshold = 1,
        mir_predicted_targets = list(mircode_symbol),
        random_seed = 1234))
})

test_that("Regression filter pass through without elastic net",{

    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        elastic.net = FALSE,
        mir_predicted_targets = list(mircode_symbol),
        random_seed = 1234)

    #check mircode_symbol how many interactions we have and intersect that with
    #the mirnas in the test set
    expect_equal(nrow(genes_miRNA_candidates$ASAP2),
                 length(intersect(names(which(mircode_symbol["ASAP2",] > 0)),
                           colnames(mir_expr))))
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A),
                 length(intersect(names(which(mircode_symbol["DNMT3A",] > 0)),
                                  colnames(mir_expr))))
})

test_that("Regression filter works with parallelization",{
    registerDoParallel(cl)

    genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
        gene_expr = gene_expr[,c("ASAP2", "DNMT3A")],
        mir_expr = mir_expr,
        mir_predicted_targets = list(mircode_symbol, targetscan_symbol),
        random_seed = 1234)

    expect_length(genes_miRNA_candidates, 2)
    expect_equal(nrow(genes_miRNA_candidates$ASAP2), 16)
    expect_equal(nrow(genes_miRNA_candidates$DNMT3A), 6)

    registerDoSEQ()
})

stopCluster(cl)