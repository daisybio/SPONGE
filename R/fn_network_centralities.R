#' Computes various node centralities
#'
#' @description Computes degree, eigenvector centrality and betweenness
#' centrality for the ceRNA interaction network induced by the results of the
#' SPONGE method
#'
#' @param sponge_result output of the sponge method
#'
#' @importFrom igraph graph.data.frame
#' @importFrom igraph eigen_centrality
#' @importFrom igraph betweenness
#' @importFrom igraph degree
#'
#' @return data.table with gene, degree, eigenvector and betweenness
#' @export
#'
#' @seealso sponge
#'
#' @examples
sponge_node_centralities <- function(sponge_result, directed = FALSE){

    network <- igraph::graph.data.frame(sponge_result)

    ev_centrality <- igraph::eigen_centrality(network,
                                              directed = directed)
    btw_centrality <- igraph::betweenness(network, directed = directed)

    page_rank <- igraph::page_rank(network, directed = directed)

    centrality_df <- data.table(gene = names(ev_centrality$vector),
                                degree = igraph::degree(network),
                                eigenvector = ev_centrality$vector,
                                betweenness = btw_centrality,
                                page_rank = page_rank$vector)
    return(centrality_df)
}

#' Computes various edge centralities
#'
#' @description Computes edge betweenness
#' centrality for the ceRNA interaction network induced by the results of the
#' SPONGE method.
#'
#' @param sponge_result
#'
#' @importFrom igraph graph.data.frame
#' @importFrom igraph E
#' @importFrom igraph edge_betweenness
#' @importFrom igraph degree
#' @importFrom tidyr separate
#'
#' @return data.table with gene, degree, eigenvector and betweenness
#' @export
#'
#' @seealso sponge
#'
#' @examples
sponge_edge_centralities <- function(sponge_result){
    directed <- FALSE

    network <- igraph::graph.data.frame(sponge_result)

    edge_labels <- attr(E(network), "vnames")
    ebtw <- edge_betweenness(network, directed = directed)
    ebtw <- data.frame(labels = edge_labels, edge_betweenness = ebtw)
    ebtw <- tidyr::separate(ebtw, col = "labels", into = c("source_gene", "target_gene"), sep = "\\|")
    return(ebtw)
}

sponge_matched_edge_betweenness <- function(cancer_network, normal_network, directed = FALSE){
    liver_cancer_edge_betweenness <- sponge_edge_betweenness(network = liver_cancer_sponge_network)
    liver_normal_edge_betweenness <- sponge_edge_betweenness(network = liver_normal_sponge_network)

    dplyr::inner_join(liver_cancer_edge_betweenness, liver_normal_edge_betweenness, by = c("source_gene", "target_gene")) %>%
        dplyr::rename(edge_betweenness_cancer = edge_betweenness.x, edge_betweenness_normal = edge_betweenness.y) %>%
        dplyr::mutate(edge_betweenness_diff = edge_betweenness_cancer - edge_betweenness_normal)
}

sponge_matched_network_centralities <- function(cancer_network, normal_network, directed = FALSE){

    cancer_centrality_df <- sponge_network_centralities(cancer_network, directed)
    normal_centrality_df <- sponge_network_centralities(normal_network, directed)

    centrality_matched <- dplyr::full_join(normal_centrality_df, cancer_centrality_df, by="gene") %>%
        dplyr::rename(eigenvector_centrality.normal = eigenvector.x,
               eigenvector_centrality.cancer = eigenvector.y,
               betweenness_centrality.normal = betweenness.x,
               betweenness_centrality.cancer = betweenness.y,
               degree.normal = degree.x,
               degree.cancer = degree.y) %>%
        mutate(eigenvector_centrality.diff = eigenvector_centrality.normal - eigenvector_centrality.cancer,
               eigenvector_rank.normal = rank(1/eigenvector_centrality.normal),
               eigenvector_rank.cancer = rank(1/eigenvector_centrality.cancer),
               eigenvector_rank.diff = eigenvector_rank.cancer - eigenvector_rank.normal,
               betweenness_centrality.diff = betweenness_centrality.normal - betweenness_centrality.cancer,
               betweenness_rank.normal = rank(1/betweenness_centrality.normal),
               betweenness_rank.cancer = rank(1/betweenness_centrality.cancer),
               betweenness_rank.diff = betweenness_rank.cancer - betweenness_rank.normal
               )

    return(centrality_matched)
}

sponge_transform_centralities_for_plotting <- function(network_centralities){
    library(dplyr)
    library(digest)
    plot.data <- bind_rows(dplyr::select(network_centralities,
                                         gene,
                                         degree = degree.cancer,
                                         eigenvector_centrality = eigenvector_centrality.cancer,
                                         betweenness_centrality = betweenness_centrality.cancer
                                         ) %>% mutate(dataset = "cancer"),
                           dplyr::select(network_centralities,
                                         gene,
                                         degree = degree.normal,
                                         eigenvector_centrality = eigenvector_centrality.normal,
                                         betweenness_centrality = betweenness_centrality.normal
                                         ) %>% mutate(dataset = "normal"))
    plot.data <- plot.data %>% mutate(color = paste("#", substr(sapply(gene, function(x) digest(x, algo = "crc32")), 1, 6), sep=""))
    return(plot.data)
}

sponge_plot_edge_centralities <- function(edge_centralities, n){
    top_x <- head(dplyr::arrange(edge_betweenness, desc(edge_betweenness)), n) %>%
        dplyr::mutate(edge = paste(source_gene, target_gene, sep="|"))
    ggplot(top_x, aes(x = reorder(edge, -edge_betweenness), y = edge_betweenness)) +
        geom_bar(stat = "identity") +
        theme_bw(base_size = base_size) +
        xlab("ceRNA interaction") +
        ylab("Edge betweenness") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

#' plot node network centralities
#'
#' @param network_centralities a result from sponge_node_centralities()
#' @param measure one of 'all', 'degree', 'ev' or 'btw'
#' @param x plot against another column in the data table, defaults to degree
#' @param top label the top x samples in the plot
#' @import ggplot2
#' @import ggrepel
#' @import digest
#' @import gridExtra
#'
#' @return a plot
#' @export
#'
#' @examples #sponge_plot_network_centralities
sponge_plot_network_centralities <- function(network_centralities,
                                             measure="all",
                                             x = "degree",
                                             top = 5,
                                             base_size = 18){

    network_centralities <- network_centralities %>% mutate(color =
        paste("#", substr(sapply(gene, function(x)
                  digest(x, algo = "crc32")), 1, 6), sep=""))

    p1 <- ggplot(network_centralities, aes(x=degree)) +
        geom_histogram(color=I("black"), fill=I("black"), alpha = 0.3)+
        theme_bw(base_size = base_size) +
        theme(strip.background = element_rect(fill="grey"))
    p2 <- ggplot(network_centralities, aes_string(x = x,
                                                  y = "eigenvector",
                                                  color= "color",
                                                  label = "gene")) +
        geom_point(alpha = 0.3) +
        ylab("eigenvector") +
        theme_bw(base_size = base_size) +
        theme(legend.position = "none") +
        geom_label_repel(data = network_centralities %>% top_n(top, eigenvector))
    p3 <- ggplot(network_centralities, aes_string(x = x,
                                                  y = "betweenness",
                                                  color = "color",
                                                  label = "gene")) +
        geom_point(alpha = 0.3) +
        ylab("betweenness") +
    theme_bw(base_size = base_size) +
        theme(legend.position = "none") +
        geom_label_repel(data = network_centralities %>% top_n(top, betweenness))
    if(measure == "degree") return(p1)
    else if(measure == "ev") return(p2)
    else if(measure == "btw") return(p3)
    else grid.arrange(p1,
                      p2 + theme(strip.background = element_blank(),
                                strip.text.x = element_blank(),
                                legend.position = "none"),
                      p3 + theme(strip.background = element_blank(),
                                  strip.text.x = element_blank(),
                                  legend.position = "none"),
                      ncol = 1)
}

sponge_plot_eigenvector_centralities_differences <- function(network_centralities, label.threshold){
    network_centralities <- network_centralities %>% mutate(color = paste("#", substr(sapply(gene, function(x) digest(x, algo = "crc32")), 1, 6), sep=""))

    ggplot(data = network_centralities,
           aes(x = eigenvector_centrality.normal,
               y = eigenvector_centrality.cancer,
               color = color, label = gene)) +
        geom_point() +
        ylab("weighted eigenvector centrality in cancer") +
        xlab("weighted eigenvector centrality in normal") +
        theme(legend.position = "none") +
        geom_label_repel(data = dplyr::filter(network_centralities, abs(eigenvector_centrality.diff) > label.threshold))
}

sponge_plot_betweenness_centralities_differences <- function(network_centralities, label.threshold){
    network_centralities <- network_centralities %>% mutate(color = paste("#", substr(sapply(gene, function(x) digest(x, algo = "crc32")), 1, 6), sep=""))

    ggplot(data = network_centralities,
           aes(x = betweenness_centrality.normal,
               y = betweenness_centrality.cancer,
               color = color, label = gene)) +
        geom_point() +
        ylab("betweenness centrality in cancer") +
        xlab("betweenness centrality in normal") +
        theme(legend.position = "none") +
        geom_label_repel(data = dplyr::filter(network_centralities, abs(betweenness_centrality.diff) > label.threshold))
}