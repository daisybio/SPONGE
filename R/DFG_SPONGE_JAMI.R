#####################################
########## PLOTS ####################
#####################################
#install.packages("ggpubr")
library(ggpubr)
library(dplyr)
library(ggplot2)

results_jami=read.table("/nfs/home/students/mhoffmann/SPONGE-TCGA/JAMI/breast_cancer_results.txt",sep="\t",header=TRUE)

results_sponge=load("/nfs/proj/SPONGE/results/breast invasive carcinoma_sponge_results_all_mirs.RData")
sponge_effects=as.data.frame(sponge_effects)

##filter only genes which both have
unique_genes=unique(unique(results_jami[["Source"]]),unique(results_jami[["Target"]]))

sponge_effects=sponge_effects[sponge_effects$geneA %in% unique_genes,]
sponge_effects=sponge_effects[sponge_effects$geneB %in% unique_genes,]

unique_genes_sponge=unique(unique(sponge_effects[["geneA"]]),unique(sponge_effects[["geneB"]]))

results_jami=results_jami[results_jami$Source %in% unique_genes_sponge,]
results_jami=results_jami[results_jami$Target %in% unique_genes_sponge,]

##filter columns which I want to have
sponge_effects_mscor=sponge_effects[, c("geneA", "geneB","mscor")]
sponge_effects_pvalue=sponge_effects[, c("geneA", "geneB","p.adj")]
sponge_effect_pcor=sponge_effects[,c("geneA","geneB","pcor")]
colnames(sponge_effects_pvalue)<-c("geneA","geneB","SPONGE_p.adj")
colnames(sponge_effect_pcor)<-c("geneA","geneB","SPONGE_pcor")

results_jami_cmi=results_jami[,c("Source","Target","CMI")]
results_jami_pvalue=results_jami[,c("Source","Target","p.adjusted")]
colnames(results_jami_cmi)<-c("geneA","geneB","CMI")
colnames(results_jami_pvalue)<-c("geneA","geneB","JAMI_p.adj")

##merge columns
data_plot_pvalue=merge(results_jami_pvalue,sponge_effects_pvalue)
data_plot_scores=merge(results_jami_cmi,sponge_effects_mscor)
data_plot_pcor_cmi=merge(results_jami_cmi,sponge_effect_pcor)

sapply(data_plot_pvalue, class)
sapply(data_plot_scores, class)
sapply(data_plot_pcor_cmi, class)

##filter significant ones

data_plot_pvalue_JAMI=data_plot_pvalue[data_plot_pvalue$JAMI_p.adj<0.05,]
unique_genes_significant_JAMI=unique(unique(data_plot_pvalue_JAMI[["geneA"]]),unique(data_plot_pvalue_JAMI[["geneB"]]))
data_plot_pvalue_SPONGE=data_plot_pvalue[data_plot_pvalue$SPONGE_p.adj<0.05,]
unique_genes_significant_SPONGE=unique(unique(data_plot_pvalue_SPONGE[["geneA"]]),unique(data_plot_pvalue_SPONGE[["geneB"]]))

data_plot_pvalue=data_plot_pvalue[data_plot_pvalue$JAMI_p.adj<0.05,]
data_plot_pvalue=data_plot_pvalue[data_plot_pvalue$SPONGE_p.adj<0.05,]
unique_genes_significant=unique(unique(data_plot_pvalue[["geneA"]]),unique(data_plot_pvalue[["geneB"]]))

##data for JAMI significant
data_plot_pcor_cmi_JAMI=data_plot_pcor_cmi[data_plot_pcor_cmi$geneA %in% unique_genes_significant_JAMI,]
data_plot_pcor_cmi_JAMI=data_plot_pcor_cmi[data_plot_pcor_cmi$geneB %in% unique_genes_significant_JAMI,]

data_plot_scores_JAMI=data_plot_scores[data_plot_scores$geneA %in% unique_genes_significant_JAMI,]
data_plot_scores_JAMI=data_plot_scores[data_plot_scores$geneB %in% unique_genes_significant_JAMI,]

##data for SPONGE significant
data_plot_pcor_cmi_SPONGE=data_plot_pcor_cmi[data_plot_pcor_cmi$geneA %in% unique_genes_significant_SPONGE,]
data_plot_pcor_cmi_SPONGE=data_plot_pcor_cmi[data_plot_pcor_cmi$geneB %in% unique_genes_significant_SPONGE,]

data_plot_scores_SPONGE=data_plot_scores[data_plot_scores$geneA %in% unique_genes_significant_SPONGE,]
data_plot_scores_SPONGE=data_plot_scores[data_plot_scores$geneB %in% unique_genes_significant_SPONGE,]

##data for both significant
data_plot_pcor_cmi=data_plot_pcor_cmi[data_plot_pcor_cmi$geneA %in% unique_genes_significant,]
data_plot_pcor_cmi=data_plot_pcor_cmi[data_plot_pcor_cmi$geneB %in% unique_genes_significant,]

data_plot_scores=data_plot_scores[data_plot_scores$geneA %in% unique_genes_significant,]
data_plot_scores=data_plot_scores[data_plot_scores$geneB %in% unique_genes_significant,]

data_plot_scores_pvalue=merge(data_plot_pvalue,data_plot_scores)
data_plot_pcor_cmi_pvalue=merge(data_plot_pvalue,data_plot_pcor_cmi)

##plot data - not for publication just to have a first look on it
plot(data_plot_pvalue$JAMI_p.adj, data_plot_pvalue$SPONGE_p.adj, main="JAMI vs SPONGE p-values",
     xlab="JAMI", ylab="SPONGE")

plot(data_plot_scores$CMI, data_plot_scores$mscor, main="CMI vs mscor",
     xlab="JAMI CMI", ylab="SPONGE mscor")

plot(data_plot_pcor_cmi$CMI, data_plot_pcor_cmi$SPONGE_pcor, main="CMI vs pcor",
     xlab="JAMI CMI", ylab="SPONGE pcor")

##plot significant ones nicely

##no pvalues
p_mscor<-ggplot(data_plot_scores_pvalue, aes(x=CMI, y=mscor)) + geom_point()+xlab("JAMI CMI") + ylab("SPONGE mscor")
p_mscor<-p_mscor+ggtitle("JAMI adj.p_value < 0.05 && SPONGE adj.p_value < 0.05")
p_mscor

p_pcor<-ggplot(data_plot_pcor_cmi_pvalue, aes(x=CMI, y=SPONGE_pcor)) + geom_point()+xlab("JAMI CMI") + ylab("SPONGE pcor")
p_pcor<-p_pcor+ggtitle("JAMI adj.p_value < 0.05 && SPONGE adj.p_value < 0.05")
p_pcor

p_mscor_JAMI<-ggplot(data_plot_scores_JAMI, aes(x=CMI, y=mscor)) + geom_point()+xlab("JAMI CMI") + ylab("SPONGE mscor")
p_mscor_JAMI<-p_mscor_JAMI+ggtitle("JAMI adj.p_value < 0.05")
p_mscor_JAMI
p_pcor_JAMI<-ggplot(data_plot_pcor_cmi_JAMI, aes(x=CMI, y=SPONGE_pcor)) + geom_point()+xlab("JAMI CMI") + ylab("SPONGE pcor")
p_pcor_JAMI<-p_pcor_JAMI+ggtitle("JAMI adj.p_value < 0.05")
p_pcor_JAMI

p_mscor_SPONGE<-ggplot(data_plot_scores_SPONGE, aes(x=CMI, y=mscor)) + geom_point()+xlab("JAMI CMI") + ylab("SPONGE mscor")
p_mscor_SPONGE<-p_mscor_SPONGE+ggtitle("SPONGE adj.p_value < 0.05")
p_mscor_SPONGE
p_pcor_SPONGE<-ggplot(data_plot_pcor_cmi_SPONGE, aes(x=CMI, y=SPONGE_pcor)) + geom_point()+xlab("JAMI CMI") + ylab("SPONGE pcor")
p_pcor_SPONGE<-p_pcor_SPONGE+ggtitle("SPONGE adj.p_value < 0.05")
p_pcor_SPONGE

##inc pvalues
p_mscor<-ggplot(data_plot_scores_pvalue, aes(x=CMI, y=mscor)) + geom_point(aes(size=JAMI_p.adj, color=SPONGE_p.adj))+xlab("JAMI CMI") + ylab("SPONGE mscor") +scale_size(trans = "reverse")
p_mscor

##check correlation
##JAMI pvalue < 0.05
#mscor
cor.test(data_plot_scores_JAMI$mscor, data_plot_scores_JAMI$CMI, method=c("pearson", "kendall", "spearman"))
#pcor
cor.test(data_plot_pcor_cmi_JAMI$SPONGE_pcor, data_plot_pcor_cmi_JAMI$CMI, method=c("pearson", "kendall", "spearman"))
##SPONGE pvalue < 0.05
#mscor
cor.test(data_plot_scores_SPONGE$mscor, data_plot_scores_SPONGE$CMI, method=c("pearson", "kendall", "spearman"))
#pcor
cor.test(data_plot_pcor_cmi_SPONGE$SPONGE_pcor, data_plot_pcor_cmi_SPONGE$CMI, method=c("pearson", "kendall", "spearman"))
##JAMI pvalue < 0.05 && SPONGE pvalue < 0.05
#mscor
cor.test(data_plot_scores$mscor, data_plot_scores$CMI, method=c("pearson", "kendall", "spearman"))
#pcor
cor.test(data_plot_pcor_cmi$SPONGE_pcor, data_plot_pcor_cmi$CMI, method=c("pearson", "kendall", "spearman"))

#####################################
#####RUN JAMI // RUN SPONGE##########
#####################################

# Title     : TODO
# Objective : TODO
# Created by: markus
# Created on: 18.05.21

#library(devtools)
#install_github("SchulzLab/RJAMI")

library(RJAMI)
library(SPONGE)
library(tidyverse)
library(rJava)

test_jvm()

########################################
##GET SPONGE BREAST CANCER CMIs ######
########################################
#load data
load(system.file("extdata","test_data.RData",package = "RJAMI"))
sapply(gene_mir_interactions_triplets,class)
sapply(gene_expr,class)
sapply(mir_expr,class)
jami_settings()
#jami_settings(tripleFormat = TRUE, considerZeros=FALSE)
#jami_settings()
result_jami <- jami(gene_miRNA_interactions = gene_mir_interactions_triplets,gene_expr = gene_expr,mir_expr = mir_expr, output_file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/result_toyset.txt")


#load interaction triplets
load("/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/results/breast invasive carcinoma_miRNA_importance_all.RData")

sponge_miRNAs=sponge_miRNAs %>% select(1:3)
colnames(sponge_miRNAs) <- c("geneA","geneB","mirnas")
sapply(sponge_miRNAs,class)
sponge_miRNAs$geneA <- as.factor(sponge_miRNAs$geneA)
sponge_miRNAs$geneB <- as.factor(sponge_miRNAs$geneB)
sponge_miRNAs$mirnas <- as.factor(sponge_miRNAs$mirnas)
sapply(sponge_miRNAs,class)
#sponge_miRNAs<-sponge_miRNAs[order(sponge_miRNAs$geneA)]

load("/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/input/breast invasive carcinoma.RData")

#delete ensg versions
gene_names <- colnames(cancer_gene_expr)
gene_names<- vapply(strsplit(gene_names,"[.]"), `[`, 1, FUN.VALUE=character(1))
colnames(cancer_gene_expr)<-gene_names

cancer_gene_expr_t <- t(cancer_gene_expr)
cancer_mir_expr_t <- t(cancer_mir_expr)

cancer_gene_expr_t<-as.data.frame(cancer_gene_expr_t)
sapply(cancer_gene_expr_t,class)
cancer_mir_expr_t<-as.data.frame(cancer_mir_expr_t)
sapply(cancer_mir_expr_t,class)

#delete genes / miRNAs which are not inside triplets
unique_triplets_geneA = unique(sponge_miRNAs[["geneA"]])
unique_triplets_geneA=as.character(unique_triplets_geneA)
unique_triplets_geneB = unique(sponge_miRNAs[["geneB"]])
unique_triplets_geneB=as.character(unique_triplets_geneB)
unique_triplets_genes=unique(c(unique_triplets_geneB, unique_triplets_geneA))

unique_mirnas=unique(sponge_miRNAs[["mirnas"]])
unique_mirnas=as.character(unique_mirnas)
expr_names=rownames(cancer_gene_expr_t)
unique_triplets_genes=as.character(unique_triplets_genes)
cancer_mir_expr_t=cancer_mir_expr_t[(row.names(cancer_mir_expr_t) %in% unique_mirnas),]
cancer_gene_expr_t=cancer_gene_expr_t[(row.names(cancer_gene_expr_t) %in% unique_triplets_genes),]

jami_settings()
jami_settings(tripleFormat = TRUE, considerZeros=FALSE)
jami_settings(tripleFormat = TRUE, considerZeros=TRUE)
#jami_settings()
result_jami <- jami(gene_miRNA_interactions = sponge_miRNAs,
                   gene_expr = cancer_gene_expr_t,
                   mir_expr = cancer_mir_expr_t, output_file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/breast_cancer_results.txt")

##AS RJAMI IS NOT WORKING AS IT SHOULD BE AFTER 2 DAYS OF TRYING I WILL USE THE ORIGINAL JAVA TOOL
write.table(cancer_gene_expr_t,file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/cancer_gene_expr_t.tsv", sep = "\t",quote=FALSE)
write.table(cancer_mir_expr_t,file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/cancer_mir_expr_t.tsv", sep = "\t",quote=FALSE)
write.table(sponge_miRNAs,file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/sponge_miRNAs.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

cancer_gene_expr_t=read.table(file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/cancer_gene_expr_t.tsv", sep="\t")
cancer_mir_expr_t=read.table(file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/cancer_mir_expr_t.tsv", sep="\t")
sponge_miRNAs=read.table(file = "/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/JAMI/sponge_miRNAs.tsv", sep="\t")
sponge_miRNAs=sponge_miRNAs[1:10000000,]

testprint("HELLO")

########################################
##GET SPONGE BREAST CANCER CMIs ######
########################################


########################################
##GET SPONGE BREAST CANCER MSCORs ######
########################################
#load data
load("/home/markus/data/SPONGEdb_TCGA/breast_cancer_analysis/input/breast invasive carcinoma.RData")

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(cancer_gene_expr[1:5,1:8])

## ---- eval=FALSE--------------------------------------------------------------
#  head(mir_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(cancer_mir_expr[1:5,1:5])

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
  mir_predicted_targets = targetscan_ensg)

## -----------------------------------------------------------------------------
genes_miRNA_candidates[1:2]

## ---- message=FALSE, warning=FALSE--------------------------------------------
ceRNA_interactions <- sponge(gene_expr = gene_expr,
                        mir_expr = mir_expr,
                        mir_interactions = genes_miRNA_candidates)

## ---- eval=FALSE--------------------------------------------------------------
#  head(ceRNA_interactions)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(ceRNA_interactions))

########################################
##GET SPONGE BREAST CANCER MSCORs ######
########################################



print("hello")

