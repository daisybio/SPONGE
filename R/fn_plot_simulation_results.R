#' Plot simulation results for different null models
#'
#' @param null_model_data the output of sponge_build_null_model
#' @import ggplot2
#' @importFrom data.table rbindlist
#' @importFrom dplyr filter
#' @return a ggplot2 object
#' @export
#'
#' @examples #null_model <- sponge_build_null_model(100, 100)
#' #sponge_plot_simulation_results(null_model)
sponge_plot_simulation_results <- function(null_model_data)
{
    test.data <- rbindlist(lapply(null_model_data, rbindlist, idcol = "k"), idcol = "m")
    test.data$k_label <- paste("cor =", test.data$k)
    test.data$m_label <- factor(paste("m =", test.data$m), levels = paste("m =", seq(1:15)))

    ggplot(dplyr::filter(test.data, m < 9)) +
        geom_density(aes(x = mscor)) +
        facet_grid(k_label ~ m_label) +
        theme_bw() +
        xlab("mscor") +
        scale_x_continuous(breaks=c(-0.15, 0, 0.15)) +
        coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(0,50))

}