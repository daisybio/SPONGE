#' Plot simulation results for different null models
#'
#' @param null_model_data the output of sponge_build_null_model
#' @importFrom data.table rbindlist
#' @return a ggplot2 object
#' @export
#'
#' @examples null_model <- sponge_build_null_model(100, 100)
#' sponge_plot_simulation_results(null_model)
sponge_plot_simulation_results <- function(null_model_data)
{
    if(!(requireNamespace("ggplot2", quietly = TRUE)))
        stop("install package ggplot2 to produce this plot")

    test.data <- rbindlist(lapply(null_model_data, rbindlist, idcol = "k"), idcol = "m")
    test.data$k_label <- paste("cor =", test.data$k)
    test.data$m_label <- factor(paste("m =", test.data$m), levels = paste("m =", seq_len(8)))

    ggplot2::ggplot(test.data) +
        ggplot2::geom_density(ggplot2::aes(x = mscor)) +
        ggplot2::facet_grid(k_label ~ m_label) +
        ggplot2::theme_bw() +
        ggplot2::xlab("mscor") +
        ggplot2::scale_x_continuous(breaks=c(-0.15, 0, 0.15)) +
        ggplot2::coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(0,50))

}