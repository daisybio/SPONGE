library(SPONGE)

#load data
load("/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/input/breast invasive carcinoma.RData")

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(cancer_gene_expr[1:5,1:8])

## ---- eval=FALSE--------------------------------------------------------------
#  head(mir_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(cancer_mir_expr[1:5,1:5])

## -- get subtypes into array
# filter NAs
cancer_gene_expr<-cancer_gene_expr[rowSums(is.na(cancer_gene_expr))==0,]
cancer_mir_expr<-cancer_mir_expr[rowSums(is.na(cancer_mir_expr))==0,]
sample_names <- row.names(cancer_gene_expr)
subtypes <- cancer_meta_data %>% dplyr::select("sampleID","_PANCAN_Cluster_Cluster_PANCAN")
subtypes <- subset(subtypes, sampleID %in% sample_names)
subtypes[is.na(subtypes)] <- "NA"
subtypes_batches <- dplyr::pull(subtypes,"_PANCAN_Cluster_Cluster_PANCAN")

## load targetscan
load("/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/input/targetscan.RData")
#filter NAs
targetscan_ensg<-targetscan_ensg[rowSums(is.na(targetscan_ensg))==0,]
targetscan_ensg<-targetscan_ensg[!(row.names(targetscan_ensg) %in% c("NA.","NA")), ]
targetscan_ensg = targetscan_ensg[-1,]

#prepare cancer_gene_expr ENSG names (delete version number)
gene_names <- colnames(cancer_gene_expr)
gene_names<- vapply(strsplit(gene_names,"[.]"), `[`, 1, FUN.VALUE=character(1))
colnames(cancer_gene_expr)<-gene_names

#regression coefficient
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
  gene_expr = cancer_gene_expr,
  mir_expr = cancer_mir_expr,
  mir_predicted_targets = targetscan_ensg,
  batches = subtypes_batches)

genes_miRNA_candidates[1:2]




