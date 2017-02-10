sponge_plot_density <- function(cancer_sponge_effects, normal_sponge_effects){
    cancer_sponge_effects$type <- "Liver cancer"
    normal_sponge_effects$type <- "Liver normal"

    p1 <- ggplot(cancer_sponge_effects) +
        facet_wrap(~type) +
        xlab("Cohen's q") +
        geom_density(aes(x = scor, fill = I("cornflowerblue")), alpha = 0.6)

    p2 <- ggplot(normal_sponge_effects) +
        facet_wrap(~type) +
        xlab("Cohen's q") +
        geom_density(aes(x = scor, fill = I("darkorchid1")), alpha = 0.6)

    p3 <- ggplot(cancer_sponge_effects) +
        xlab("p-value") +
        geom_density(aes(x = scor_p, fill = I("cornflowerblue")), alpha = 0.6)

    p4 <- ggplot(normal_sponge_effects) +
        xlab("p-value") +
        geom_density(aes(x = scor_p, fill = I("darkorchid1")), alpha = 0.6)

    p5 <- ggplot(cancer_sponge_effects) +
        xlab("adjusted p-value (BH)") +
        geom_density(aes(x = scor_p.adj, fill = I("cornflowerblue")), alpha = 0.6)

    p6 <- ggplot(normal_sponge_effects) +
        xlab("ajusted p-value (BH)") +
        geom_density(aes(x = scor_p.adj, fill = I("darkorchid1")), alpha = 0.6)

    grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
}

sponge_plot_heatmap <- function(data, interactive=T, show = "scor"){
    if(require(d3heatmap) && interactive){
        sponge.matrix <- xtabs(p.adj ~ source_gene + target_gene, data = data)
        d3heatmap(sponge.matrix, dendrogram = "none", symm=T)
    }
    else if(require(ggplot2))
    {
        if(interactive) warn("Library d3heatmap not found but required for interactive plotting")
        data$significance <- ""
        data[data$p.adj < 0.05, "significance"] <- "*"
        data[data$p.adj < 0.01, "significance"] <- "**"
        data[data$p.adj < 0.001, "significance"] <- "***"
        ggplot(data = data, aes_string(fill = show, x = "source_gene", y = "target_gene")) +
            geom_tile() +
            theme_bw() +
            geom_text(label=data$significance)
    }
    else{
        stop("No suitable plotting library found. Please install d3heatmap for interactive and ggplot2 for static heatmaps")
    }
}

sponge_plot_boxplot <- function(data){
    if(!require(ggplot2)) stop("library ggplot2 needs to be installed for this plot")

    ggplot(data = data, aes(x = source_gene, y = scor)) + geom_boxplot(fill = "skyblue", aes(outlier.color = p.adj)) + theme_bw()
}

sponge_network <- function(sponge.data,
                           mir.data,
                           target.genes = NULL,
                           cerna.p.val.threshold = 0.05,
                           show.sponge.interaction = TRUE,
                           show.mirnas = c("none", "all", "shared"),
                           replace.mirna.with.name = TRUE,
                           min.interactions = 3){
    library(foreach)
    library(iterators)
    library(dplyr)
    sponge.data <- filter(sponge.data, p.val < cerna.p.val.threshold)

    #genes <- unique(c(as.character(sponge.data$source_gene), as.character(sponge.data$target_gene)))
    genes <- unique(c(as.character(sponge.data$source_gene), as.character(sponge.data$target_gene)))

    edges <- NULL
    nodes <- NULL

    #sponge genes edges
    sponge.edges <- with(sponge.data, {
        result <- data.frame(from=source_gene, to=target_gene, width = scor, color="red")

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
            gene_mirnas <- mir.data[[as.character(gene)]]
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

sponge_plot_network <- function(nodes,
                                edges,
                                layout="layout.fruchterman.reingold",
                                force.directed = FALSE){
    library(visNetwork)

    if(nrow(edges) < 5000){
        plot <- visNetwork(nodes, edges)
        plot <- plot %>% visIgraphLayout(layout = layout, type="full", physics = force.directed)
        plot <- plot %>% visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
        plot <- plot %>% visNodes(font = list(size = 32))
        plot <- plot %>% visEdges(color = list(opacity = 1))

        return(plot)
    }
    else{
        warning("With more than 5000 edges, this plot is omitted for performance reasons")
    }
    return(edges)
}