library(ggplot2)
library(reshape2)
# plot expression of miRNAs in module for specific RNA (no negative values!)
plot_miRNAs_per_gene <- function(target, genes_miRNA_candidates, mi_rna_expr, meta, log_transform = T, pseudocount = 1e-3, unit = "counts") {
  # extract miRNAs associated with given RNA
  candidates <- genes_miRNA_candidates[[target]]
  miRNAs <- candidates[["mirna"]]
  miRNA_expr <- mi_rna_expr[,miRNAs]
  
  # sum all expressions for the given conditions
  miRNA_expr <- merge(miRNA_expr, meta[,c("disease_id", "disease_type")], by.x = 0, by.y = "disease_id")
  miRNA_expr <- data.frame(miRNA_expr, row.names = 1)
  miRNA_expr <- aggregate(miRNA_expr[,1:(ncol(miRNA_expr)-1)], list(condition=miRNA_expr$disease_type), FUN=sum)
  miRNA_expr <- t(data.frame(miRNA_expr, row.names = 1))
  miRNA_expr <- melt(miRNA_expr, id.vars = "x", varnames = c("miRNA", "variable"))
  miRNA_expr$miRNA <- gsub("\\.", "-", miRNA_expr$miRNA)
  
  # label
  label <- paste0("miRNA ", unit, " over samples")
  # log transform if given
  if (log_transform){
    miRNA_expr$value <- log10(miRNA_expr$value + pseudocount)
    label <- paste0("log10 + ", pseudocount, " ", label)
  }
  # two bar plots in one
  bars <- ggplot(miRNA_expr, aes(y=miRNA, x=value, fill=variable)) +
    ggtitle(target) +
    geom_bar(stat='identity', position='dodge') +
    xlab(label) + 
    scale_fill_discrete(name = "Condition")
  return(bars)
}

library(pheatmap)
# TG heat map
plot_target_gene_expressions <- function(target, target_genes, gene_expression, meta, log_transform = T, pseudocount = 1, gtf_raw = NA, annotation = NA) {
  # get target expression
  target_expression <- gene_expression[,target, drop = F]
  # get target genes expressions
  target_genes_expression <- gene_expression[,target_genes]
  # combine into one
  data <- t(cbind(target_expression, target_genes_expression))
  # save all conditions of samples
  conditions <- unique(meta$disease_type)
  
  # convert to hgnc if gtf is given
  if (!is.na(gtf_raw)) {
    print("reading GTF file and converting to HGNC symbols...")
    gtf <- rtracklayer::readGFF(gtf)
    gene.ens.all <- unique(gtf[!is.na(gtf$transcript_id),c("gene_id", "gene_name")])
    colnames(gene.ens.all) <- c("ensembl_gene_id", "hgnc_symbol")
    rownames(gene.ens.all) <- gene.ens.all$ensembl_gene_id
    # convert
    annotation = gene.ens.all
  }
  if (!is.na(annotation)){
    ensgs <- rownames(data)
    ensgs <- merge(ensgs, annotation, by = 1, all.x = T)
    ensgs[!is.na(ensgs$hgnc_symbol),"x"] <- ensgs[!is.na(ensgs$hgnc_symbol),"hgnc_symbol"]
    rownames(data) <- ensgs$x
  }
  # heat color scheme
  colors <- c(colorRampPalette(c("blue", "orange"))(100), colorRampPalette(c("orange", "red"))(100))
  # annotation colors
  annotation.colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # annotation.colors <- met.brewer(palette, n = length(conditions))
  # name colors
  names(annotation.colors) <- conditions
  # create annotation for heat map
  df <- data.frame(meta[,c("disease_id", "disease_type")], row.names = 1)
  data <- data+1
  pheatmap::pheatmap(data, treeheight_row = 0, treeheight_col = 0,
                    show_colnames = F, cluster_rows = T, cluster_cols = T,
                   color = colors, annotation_col = df,
                  annotation_colors = list(condition=annotation.colors), main = target, fontsize_row = 10)
}