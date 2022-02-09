##CONTENT
## 1. SPONGE RUN FOR circRNA-miRNA-mRNA network
## 2. SPONGE RUN FOR circRNA-miRNA network

#################################################
## 1. SPONGE RUN FOR circRNA-miRNA-mRNA network##
#################################################

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(SPONGE)
library(dplyr)
library(reshape2)
library(doParallel)
registerDoParallel(cores=4)
#a. load circRNA expression
circRNA_expression<-read.csv("/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/input_files_SPONGE/circRNA_counts_filtered_transformed.tsv", sep='\t')
#b. load mirNA expression
#official
#mirna_expression<-read.csv("/home/markus/data/SPONGEdb_TCGA/circRNA/input_files_SPONGE/TGCT_NEW_cohort_n114.sum_read_counts_per_MIMAT.NAs_to_zeros.txt",sep='\t')
#normalized counts octavias pipeline
mirna_expression<-read.csv("/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/input_files_SPONGE/miRNA_counts_filtered.tsv",sep='\t', row.names = 1)
#c. load binding site prediction
circRNA_binding_site_prediction<-read.csv("/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/input_files_SPONGE/bindsites_25%_filtered_sponge_transformed.tsv", sep='\t')
#d. check if mirna and circRNA have the same samples, all others delete + order + check again afterwards + transform matrix +convert to matrix
col_names_circRNA <- colnames(circRNA_expression)
col_names_miRNA <- colnames(mirna_expression)
intersect(col_names_circRNA, col_names_miRNA)
mirna_expression<-mirna_expression %>% select(one_of(c(col_names_circRNA)))
circRNA_expression<-circRNA_expression %>% select(one_of(c(col_names_miRNA)))
mirna_expression <- mirna_expression[,order(colnames(mirna_expression))]
circRNA_expression <- circRNA_expression[,order(colnames(circRNA_expression))]
col_names_circRNA <- colnames(circRNA_expression)
col_names_miRNA <- colnames(mirna_expression)
circRNA_expression<-t(circRNA_expression)
mirna_expression<-t(mirna_expression)
#circRNA_expression<-as.matrix(circRNA_expression)
#mirna_expression<-as.matrix(mirna_expression)
circRNA_binding_site_prediction<-as.matrix(circRNA_binding_site_prediction)
colnames(circRNA_binding_site_prediction)<-gsub('\\.','-',colnames(circRNA_binding_site_prediction))

intersect(colnames(mirna_expression[1:10,1:10]), colnames(circRNA_binding_site_prediction))

#e. predict circRNA-miRNA candidates
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
  gene_expr = circRNA_expression,
  mir_expr = mirna_expression,
  mir_predicted_targets = circRNA_binding_site_prediction,
  coefficient.threshold=0)

circRNA_miRNA_candidates<-genes_miRNA_candidates

save(circRNA_miRNA_candidates, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_01_circRNA_miRNA_candidates.RDATA")

#f. predict ceRNA_interactions
ceRNA_interactions <- sponge(gene_expr = circRNA_expression,
                             mir_expr = mirna_expression,
                             mir_interactions = circRNA_miRNA_candidates)

save(ceRNA_interactions, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_02_ceRNA_interactions.RDATA")


#g. predict null model
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, number_of_samples = nrow(gene_expr))

save(mscor_null_model, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_03_mscor_null_model.RDATA")

#h. plot simulation results
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_04_circRNA_miRNA_sponge_plot_simulation_results.png")
sponge_plot_simulation_results(mscor_null_model)
#i. calculcate sign. ceRNA interactions
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions,
                                                   null_model = mscor_null_model)

save(ceRNA_interactions_sign, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_05_ceRNA_interactions_sign.RDATA")

#j. calculate ceRNA interactions fdr
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < 0.2),]
save(ceRNA_interactions_fdr, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_06_ceRNA_interactions_fdr.RDATA")

#k. plot network
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_07_circRNA_miRNA_sponge_plot_network.png")
sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)

#l. calculate network centralities
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
save(network_centralities, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_08_network_centralities.RDATA")

#m. calculate weighted network centralities
ceRNA_interactions_fdr_weight <- ceRNA_interactions_fdr
ceRNA_interactions_fdr_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
save(weighted_network_centralities, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_09_weighted_network_centralities.RDATA")

#n. plot network centralities
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_10_circRNA_miRNA_sponge_plot_network_centralities.png")
sponge_plot_network_centralities(weighted_network_centralities, measure = "btw", top = 1)

#o. more covariance matrix
more_covariance_matrices <- sample_zero_mscor_cov(m = 1,
                                                  number_of_solutions = 10,
                                                  gene_gene_correlation = 0.5)
save(more_covariance_matrices, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_11_more_covariance_matrices.RDATA")

#p. mscor_coefficients

mscor_coefficients <- sample_zero_mscor_data(cov_matrices = more_covariance_matrices,
                                             number_of_samples = 200, number_of_datasets = 100)
save(mscor_coefficients, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_12_mscor_coefficients.RDATA")

#q. plot mscor coefficients
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_10_circRNA_miRNA_sponge_plot_mscor_coefficients.png")
hist(unlist(mscor_coefficients), main = "Random distribution of mscor coefficients")

##############################################
## 2. SPONGE RUN FOR circRNA-miRNA network####
##############################################

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(SPONGE)
library(dplyr)
library(reshape2)
library(doParallel)
registerDoParallel(cores=4)
#a. load circRNA expression
circRNA_expression<-read.csv("/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/input_files_SPONGE/circRNA_counts_filtered_transformed.tsv", sep='\t')
#b. load mirNA expression
#official
#mirna_expression<-read.csv("/home/markus/data/SPONGEdb_TCGA/circRNA/input_files_SPONGE/TGCT_NEW_cohort_n114.sum_read_counts_per_MIMAT.NAs_to_zeros.txt",sep='\t')
#normalized counts octavias pipeline
mirna_expression<-read.csv("/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/input_files_SPONGE/miRNA_counts_filtered.tsv",sep='\t', row.names = 1)
#c. load binding site prediction
circRNA_binding_site_prediction<-read.csv("/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/input_files_SPONGE/bindsites_25%_filtered_sponge_transformed.tsv", sep='\t')
#d. check if mirna and circRNA have the same samples, all others delete + order + check again afterwards + transform matrix +convert to matrix
col_names_circRNA <- colnames(circRNA_expression)
col_names_miRNA <- colnames(mirna_expression)
intersect(col_names_circRNA, col_names_miRNA)
mirna_expression<-mirna_expression %>% select(one_of(c(col_names_circRNA)))
circRNA_expression<-circRNA_expression %>% select(one_of(c(col_names_miRNA)))
mirna_expression <- mirna_expression[,order(colnames(mirna_expression))]
circRNA_expression <- circRNA_expression[,order(colnames(circRNA_expression))]
col_names_circRNA <- colnames(circRNA_expression)
col_names_miRNA <- colnames(mirna_expression)
circRNA_expression<-t(circRNA_expression)
mirna_expression<-t(mirna_expression)
#circRNA_expression<-as.matrix(circRNA_expression)
#mirna_expression<-as.matrix(mirna_expression)
circRNA_binding_site_prediction<-as.matrix(circRNA_binding_site_prediction)
colnames(circRNA_binding_site_prediction)<-gsub('\\.','-',colnames(circRNA_binding_site_prediction))

intersect(colnames(mirna_expression[1:10,1:10]), colnames(circRNA_binding_site_prediction))

#e. predict circRNA-miRNA candidates
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
  gene_expr = circRNA_expression,
  mir_expr = mirna_expression,
  mir_predicted_targets = circRNA_binding_site_prediction,
  coefficient.threshold=0)

circRNA_miRNA_candidates<-genes_miRNA_candidates

save(circRNA_miRNA_candidates, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_01_circRNA_miRNA_candidates.RDATA")

#f. predict ceRNA_interactions
ceRNA_interactions <- sponge(gene_expr = circRNA_expression,
                             mir_expr = mirna_expression,
                             mir_interactions = circRNA_miRNA_candidates)

save(ceRNA_interactions, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_02_ceRNA_interactions.RDATA")


#g. predict null model
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, number_of_samples = nrow(gene_expr))

save(mscor_null_model, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_03_mscor_null_model.RDATA")

#h. plot simulation results
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_04_circRNA_miRNA_sponge_plot_simulation_results.png")
sponge_plot_simulation_results(mscor_null_model)
#i. calculcate sign. ceRNA interactions
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions,
                                                   null_model = mscor_null_model)

save(ceRNA_interactions_sign, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_05_ceRNA_interactions_sign.RDATA")

#j. calculate ceRNA interactions fdr
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < 0.2),]
save(ceRNA_interactions_fdr, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_06_ceRNA_interactions_fdr.RDATA")

#k. plot network
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_07_circRNA_miRNA_sponge_plot_network.png")
sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)

#l. calculate network centralities
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
save(network_centralities, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_08_network_centralities.RDATA")

#m. calculate weighted network centralities
ceRNA_interactions_fdr_weight <- ceRNA_interactions_fdr
ceRNA_interactions_fdr_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
save(weighted_network_centralities, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_09_weighted_network_centralities.RDATA")

#n. plot network centralities
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_10_circRNA_miRNA_sponge_plot_network_centralities.png")
sponge_plot_network_centralities(weighted_network_centralities, measure = "btw", top = 1)

#o. more covariance matrix
more_covariance_matrices <- sample_zero_mscor_cov(m = 1,
                                                  number_of_solutions = 10,
                                                  gene_gene_correlation = 0.5)
save(more_covariance_matrices, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_11_more_covariance_matrices.RDATA")

#p. mscor_coefficients

mscor_coefficients <- sample_zero_mscor_data(cov_matrices = more_covariance_matrices,
                                             number_of_samples = 200, number_of_datasets = 100)
save(mscor_coefficients, file = "/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_12_mscor_coefficients.RDATA")

#q. plot mscor coefficients
png(file="/nfs/home/students/mhoffmann/SPONGE-TCGA/sponge_circRNA/circRNA_SPONGE_RUN/results/circRNA_miRNA/step_10_circRNA_miRNA_sponge_plot_mscor_coefficients.png")
hist(unlist(mscor_coefficients), main = "Random distribution of mscor coefficients")