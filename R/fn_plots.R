#' Prepare a sponge network for plotting
#' @import foreach
#' @import iterators
#' @importFrom dplyr filter
#' @param sponge_result ceRNA interactions as produced by the sponge method.
#' @param mir_data miRNA interactions as produced by sponge_gene_miRNA_interaction_filter
#' @param target.genes a character vector  to select a subset of genes
#' @param show.sponge.interaction whether to connect ceRNAs
#' @param show.mirnas one of none, shared, all
#' @param replace.mirna.with.name uses mirbase to replace mimats with names
#' @param min.interactions minimum degree of a gene to be shown
#'
#' @return a list of nodes and edges
#' @export
#'
#' @examples sponge_network(ceRNA_interactions, mir_interactions)
sponge_network <- function(sponge_result,
                           mir_data,
                           target.genes = NULL,
                           show.sponge.interaction = TRUE,
                           show.mirnas = "none",
                           replace.mirna.with.name = TRUE,
                           min.interactions = 3){

    genes <- unique(c(as.character(sponge_result$geneA), as.character(sponge_result$geneB)))

    edges <- NULL
    nodes <- NULL

    #sponge genes edges
    sponge.edges <- with(sponge_result, {
        result <- data.frame(from=geneA, to=geneB, width = mscor, color="red")

        result[which(result$to %in% genes | result$from %in% genes),]
    })

    if(show.sponge.interaction) edges <- rbind(edges, sponge.edges)

    #sponge gene nodes
    sponge.nodes <- data.frame(id = genes, label = genes, shape = "square", color="darkgreen", type=FALSE, value=1)
    nodes <- rbind(nodes, sponge.nodes)

    if(!is.null(target.genes)){
        nodes[which(nodes$id %in% target.genes), "shape"] <- "square"
        nodes[which(nodes$id %in% target.genes), "color"] <- "blue"
    }

    #mirna nodes and edges
    if(show.mirnas != "none"){

        mirna.edges <- foreach(gene = nodes$id, .combine=rbind) %do% {
            gene_mirnas <- mir_data[[as.character(gene)]]
            gene_mirnas <- dplyr::filter(gene_mirnas, coefficient < 0)
            foreach(mir = iter(gene_mirnas, by="row"), .combine = rbind) %do% {
                with(mir, { data.frame(from = gene, to = mirna, width = abs(log2(coefficient)), color="blue")})
            }
        }

        if(show.mirnas == "shared")
        {
            #consider only miRNAs shared with a source gene
            mirnas <- mirna.edges[which(mirna.edges$from %in% genes),"to"]

            #count mirnas that appear form at least x edges
            mirnas <- intersect(mirnas, names(which(table(as.character(mirna.edges$to)) > min.interactions)))

            mirna.edges <- mirna.edges[which(mirna.edges$to %in% mirnas),]
        }

        else mirnas <- unique(as.character(mirna.edges$to))

        if(length(mirnas) > 0){
            nodes <- rbind(nodes, data.frame(id = mirnas, label = mirnas, shape = "triangle", color="darkblue", type=TRUE, value=1))
            edges <- rbind(edges, mirna.edges)

            #replace mirna with name
            if(replace.mirna.with.name){

                nodes$label <- as.character(nodes$label)
                nodes[which(nodes$type), "label"] <- as.character(fn_map_mimat_to_mir(as.character(nodes[which(nodes$type), "id"])))
            }
        }
        else{
            stop("No miRNAs found that match all criteria")
        }
    }

    #filter out orphan nodes
    nodes <- nodes[which(nodes$id %in% c(as.character(edges$to), as.character(edges$from))),]

    return(list(nodes=nodes, edges=edges))
}

#' Plot a sponge network
#'
#' @param sponge_result ceRNA interactions as produced by the sponge method.
#' @param mir_data miRNA interactions as produced by sponge_gene_miRNA_interaction_filter
#' @param layout one of the layout methods supported in the visNetwork package
#' @param force.directed whether to produce a force directed network, gets
#' slow for large networks
#' @param ... further params for sponge_network
#' @import visNetwork
#'
#' @return shows a plot
#' @export
#'
#' @examples sponge_plot_network(ceRNA_interactions, mir_interactions)
sponge_plot_network <- function(sponge_result, mir_data,
                                layout="layout.fruchterman.reingold",
                                force.directed = FALSE, ...){
    network <- sponge_network(sponge_result, mir_data, ...)
    nodes <- network$nodes
    edges <- network$edges

    if(nrow(edges) < 10000){
        plot <- visNetwork(nodes, edges)
        plot <- plot %>% visIgraphLayout(layout = layout, type="full", physics = force.directed)
        plot <- plot %>% visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
        plot <- plot %>% visNodes(font = list(size = 32))
        plot <- plot %>% visEdges(color = list(opacity = 1))

        return(plot)
    }
    else{
        warning("With more than 10000 edges, this plot is omitted for performance reasons")
    }
    return(edges)
}