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
sponge_node_centralities <- function(sponge_result){
    directed <- FALSE

    network <- igraph::graph.data.frame(sponge_result)

    ev_centrality <- igraph::eigen_centrality(network, directed = directed)
    btw_centrality <- igraph::betweenness(network, directed = directed)
    centrality_df <- data.table(gene = names(ev_centrality$vector),
                                degree = igraph::degree(network),
                                eigenvector = ev_centrality$vector,
                                betweenness = btw_centrality)
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
        theme_bw() +
        xlab("ceRNA interaction") +
        ylab("Edge betweenness") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

sponge_plot_network_centralities <- function(network_centralities, measure="all"){
    library(ggplot2)
    library(ggrepel)

    plot.data <- sponge_transform_centralities_for_plotting(network_centralities)

    p1 <- ggplot(plot.data, aes(x=degree)) +
        geom_histogram(color=I("black"), fill=I("black"), alpha = 0.3)+
        facet_wrap(~dataset, scales = "free_x") +
        xlab("node degree") +
        theme(strip.background = element_rect(fill="grey"))
    p2 <- ggplot(plot.data, aes(x = degree, y = eigenvector_centrality, color= color, label = gene)) +
        geom_point(alpha = 0.3) +
        facet_wrap(~dataset, scales = "free_x") +
        xlab("node degree") +
        ylab("weighted eigenvector centrality") +
        theme(legend.position = "none") +
        geom_label_repel(data = dplyr::group_by(plot.data, dataset) %>% top_n(5, eigenvector_centrality))
    p3 <- ggplot(plot.data, aes(x = degree, y = betweenness_centrality, color = color, label = gene)) +
        geom_point(alpha = 0.3) +
        facet_wrap(~dataset, scales = "free_x") +
        xlab("node degree") +
        ylab("betweenness centrality") +
        theme(legend.position = "none") +
        geom_label_repel(data = dplyr::group_by(plot.data, dataset) %>% top_n(5, betweenness_centrality))
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

sponge_plot_top_centralities <- function(network_centralities, top = 50,
                                         known.sponge.genes = c("ESR1", "CD44", "LIN28B", "HULC", "KRAS1P", "HSUR1", "HSUR2", "BRAFP1", "VCAN", "LINCMD1", "H19"),
                                         known.cancer.genes = c("TP53", "ESR1", "CD44", "KRAS"),
                                         only=""){

    plot.bars <- function(plot.data, value, ylabel){
        ggplot(plot.data) +
            geom_bar(aes_string(x = "gene", y = value, fill = "ceRNA", color = "cancer"),
                     stat="identity",
                     alpha = 0.6) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_manual(values = c("black","red")) +
            scale_color_manual(values = c("black","green"),
                               guide = guide_legend(title = "cancer gene")) +
            facet_wrap(~plot)
    }

    network_centralities$ceRNA = "novel"
    network_centralities[which(network_centralities$gene %in% known.sponge.genes), "ceRNA"] <- "known"
    network_centralities$ceRNA <- factor(network_centralities$ceRNA, levels = c("novel", "known"))

    network_centralities$cancer = FALSE
    network_centralities[which(network_centralities$gene %in% known.cancer.genes), "cancer"] <- TRUE
    network_centralities$cancer <- factor(network_centralities$cancer, levels = c(FALSE, TRUE))

    top_eigenvector_centrality.cancer <- head(dplyr::arrange(network_centralities, desc(eigenvector_centrality.cancer)), top)
    top_eigenvector_centrality.cancer <- within(top_eigenvector_centrality.cancer, gene <- factor(gene, levels=unique(as.character(gene))))
    top_eigenvector_centrality.cancer$plot <- "Liver cancer"

    top_eigenvector_centrality.normal <- head(dplyr::arrange(network_centralities, desc(eigenvector_centrality.normal)), top)
    top_eigenvector_centrality.normal <- within(top_eigenvector_centrality.normal, gene <- factor(gene, levels=unique(as.character(gene))))
    top_eigenvector_centrality.normal$plot <- "Liver normal"

    top_betweenness_centrality.cancer <- head(dplyr::arrange(network_centralities, desc(betweenness_centrality.cancer)), top)
    top_betweenness_centrality.cancer <- within(top_betweenness_centrality.cancer, gene <- factor(gene, levels=unique(as.character(gene))))
    top_betweenness_centrality.cancer$plot <- "Liver cancer"

    top_betweenness_centrality.normal <- head(dplyr::arrange(network_centralities, desc(betweenness_centrality.normal)), top)
    top_betweenness_centrality.normal <- within(top_betweenness_centrality.normal, gene <- factor(gene, levels=unique(as.character(gene))))
    top_betweenness_centrality.normal$plot <- "Liver normal"

    p1 <- plot.bars(top_eigenvector_centrality.cancer,
                    value = "eigenvector_centrality.cancer") +
          ylab("weighted eigenvector centrality")

    p2 <- plot.bars(top_eigenvector_centrality.normal,
                    value = "eigenvector_centrality.normal") +
          ylab("weighted eigenvector centrality")

    p3 <- plot.bars(top_betweenness_centrality.cancer,
                    value = "betweenness_centrality.cancer") +
          ylab("betweenness centrality")

    p4 <- plot.bars(top_betweenness_centrality.normal,
                    value = "betweenness_centrality.normal") +
          ylab("betweenness centrality")

    if(only == "cancer")
        grid.arrange(p1, p3)
    else if(only == "normal")
        grid.arrange(p2, p4)
    else if(only == "eigenvector"){
        grid.arrange(p1, p2)
    }
    else if(only == "betweenness"){
        grid.arrange(p3, p4)
    }
    else grid.arrange(p1, p2, p3, p4)
}