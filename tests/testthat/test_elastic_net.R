library(SPONGE)

context("TEST elastic net related functions")

test_that("TEST computing residual sum of squares", {
    set.seed(1234)
    model <- cv.glmnet(mir_expr, gene_expr[,2], alpha = 0.5)

    expect_equal(fn_get_rss(model, mir_expr, gene_expr[,2]), 7.493247,
                 tolerance = 1e-6)
})


test_that("TEST extract model coefficients", {
    set.seed(1234)
    model <- cv.glmnet(mir_expr, gene_expr[,2], alpha = 0.5)

    expect_equal(nrow(fn_get_model_coef(model)), 49)
    expect_equal(ncol(fn_get_model_coef(model)), 2)
    expect_equal(mean(fn_get_model_coef(model)$coefficient), -0.01653435)
})

test_that("TEST elastic net", {
    set.seed(1234)
    result <- fn_elasticnet(mir_expr, gene_expr[,2], alpha.step = 0.5)

    expect_equal(attr(result, "class"), "cv.glmnet")
    expect_equal(result$lambda.min, 0.03065094, tolerance = 1e-7)
    expect_equal(result$lambda[1], 1.153988, tolerance = 1e-6)
    expect_length(result$lambda, 71)
})

test_that("TEST F test", {
    set.seed(1234)
    model <- fn_elasticnet(mir_expr, gene_expr[,2], alpha.step = 0.5)
    result <- fn_gene_miRNA_F_test(model = model, g_expr = gene_expr[,2],
                                   m_expr = mir_expr)
    expect_equal(nrow(result), 44)
    expect_equal(ncol(result), 4)
    expect_equal(mean(result$p.adj), 0.1525966, tolerance = 1e-7)
})