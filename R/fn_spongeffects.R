#' Calculate spongEffects modules and immunedeconvolution scores for each cancer type
#'
#' @importFrom tidyverse
#' @importFrom caret
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom dplyr
#' @source ./fn_spongeffects_utility.R
#'
#' @param ceRNA_expression_data ceRNA expression data (same structure as input
#' for SPONGE)
#'
#'
#' @param cores parallel cores used in doParallel package (default = 1)
#' @return A vector with shared RNAs of the two genes.

fn_spongeffects_doTraining <- function(ceRNA_expression_data,
                                       cores = 1){
  # Define cores for parallel computing
  no_cores <- cores
  cl <- makePSOCKcluster(25)
  registerDoParallel(cl)

}

# Tidy and prepare input dataset
TCGA.expr.tumor <- cancer_gene_expr %>% t() %>% as.data.frame()
TCGA.expr.normal <-normal_gene_expr %>% t() %>% as.data.frame()
TCGA.meta.normal <- normal_meta_data %>%
  dplyr::filter(sampleID %in% colnames(TCGA.expr.normal))
TCGA.meta.tumor <- cancer_meta_data %>%
  dplyr::filter(sampleID %in% colnames(TCGA.expr.tumor))

TCGA.meta.tumor$PATIENT_ID <- substr(TCGA.meta.tumor$sampleID,1,nchar(TCGA.meta.tumor$sampleID)-3)

# Load metadata with grading information and keep only samples with StageI, II, III, and IV
MetaData_TCGA_BRCA_Complete <- read.delim("./Data/Metadata/BRCA_data_clinical_patient.txt", comment.char="#") %>%
  dplyr::filter(!(SUBTYPE %in% c("")) & PATIENT_ID %in% TCGA.meta.tumor$PATIENT_ID) %>%
  dplyr::filter(AJCC_PATHOLOGIC_TUMOR_STAGE %in% c('STAGE I', 'STAGE IA', 'STAGE IB', 'STAGE II', 'STAGE IIA', 'STAGE IIB', 'STAGE III', 'STAGE IIIA', 'STAGE IIIB', 'STAGE IIIC', 'STAGE IV')) %>%
  mutate(Grading = ifelse(AJCC_PATHOLOGIC_TUMOR_STAGE %in% c('STAGE I', 'STAGE IA', 'STAGE IB'), "Stage_I",
                          ifelse(AJCC_PATHOLOGIC_TUMOR_STAGE %in% c('STAGE II', 'STAGE IIA', 'STAGE IIB'), "Stage_II",
                                 ifelse(AJCC_PATHOLOGIC_TUMOR_STAGE %in% "STAGE IV", "Stage_IV", "Stage_III"))))

TCGA.meta.tumor <- TCGA.meta.tumor[match(MetaData_TCGA_BRCA_Complete$PATIENT_ID, TCGA.meta.tumor$PATIENT_ID), ] %>%
  left_join(MetaData_TCGA_BRCA_Complete, .by = PATIENT_ID)
TCGA.meta.tumor$SUBTYPE <- gsub("BRCA_","", TCGA.meta.tumor$SUBTYPE)
TCGA.meta.tumor$SUBTYPE <- factor(TCGA.meta.tumor$SUBTYPE, levels = c("LumA",  "LumB",  "Her2",  "Basal", "Normal"))

TCGA.expr.normal <- TCGA.expr.normal[, match(TCGA.meta.normal$sampleID, colnames(TCGA.expr.normal))]
TCGA.expr.tumor <- TCGA.expr.tumor[, match(TCGA.meta.tumor$sampleID, colnames(TCGA.expr.tumor))]

rownames(TCGA.expr.tumor) <- gsub("\\..*","",rownames(TCGA.expr.tumor))
rownames(TCGA.expr.normal) <- gsub("\\..*","",rownames(TCGA.expr.normal))


# Load METABRIC Data ----------------------------------------------------------
# Cohort 1 and 2

METABRIC.expr <- read.delim("./Data/Validation/brca_metabric/data_expression_median.txt",
                            header=T) %>% dplyr::select(-Entrez_Gene_Id)
METABRIC.expr <- METABRIC.expr[Biobase::isUnique(METABRIC.expr$Hugo_Symbol), ]

METABRIC.expr <- Convert_Gene_Names(METABRIC.expr)
colnames(METABRIC.expr) <- gsub("\\.", "-", colnames(METABRIC.expr))

METABRIC.meta <- read.delim("./Data/Validation/brca_metabric/data_clinical_patient.txt", comment.char="#") %>%
  filter(CLAUDIN_SUBTYPE %in% c("LumA", "LumB", "Basal", "Her2", "Normal") & PATIENT_ID %in% colnames(METABRIC.expr))
METABRIC.meta$CLAUDIN_SUBTYPE <- factor(METABRIC.meta$CLAUDIN_SUBTYPE, levels = c("LumA",  "LumB",  "Her2",  "Basal", "Normal"))

METABRIC.expr <- METABRIC.expr[, match(METABRIC.meta$PATIENT_ID, colnames(METABRIC.expr))]

# Define network ----------------------------------------------------------
# List Sponge networks
SpogeNetwork.files <- list.files("./Data/Sponge_Networks/", pattern= ".RData")
SpogeNetwork.files <- gsub("_sponge_results.RData", "", SpogeNetwork.files)
SpogeNetwork.files <- SpogeNetwork.files[ SpogeNetwork.files != "HepG2"]

cancerType <- "breast invasive carcinoma"
load(paste0("./Data/Sponge_Networks/", cancerType, "_sponge_results.RData"))

#Filter SPONGE network for significant edges
Sponge.filtered <- sponge_effects %>%
  filter_network(mscor.treshold =  .1, padj.treshold = .01)

# Calculate weighted centrality scores and add them to the ones present in SpongeDB
cancerType.Centrality <- gsub(" ", "_", cancerType)
Node.Centrality <- read.csv(paste0("./Data/NetworkAnalysis/", cancerType.Centrality, "_networkAnalysis.csv"), sep = " ")

Nodes <- Weighted_degree(Sponge.filtered, undirected = T, Alpha = 1)

Node.Centrality <- Node.Centrality %>%
  mutate(Weighted_Degree = Nodes$Weighted_degree[match(Node.Centrality$gene, Nodes$Nodes)])

# Filter for top1000 most central genes
ensembl=useMart("ensembl")
ensemblH = useDataset("hsapiens_gene_ensembl",mart=ensembl)
Human.lncRNA = getBM(attributes=c("ensembl_gene_id","gene_biotype","description"), mart=ensemblH, useCache = FALSE) %>%
  dplyr::filter(gene_biotype == "lncRNA")

Node.Centrality.weighted <- Node.Centrality %>%
  dplyr::filter(gene %in% Human.lncRNA$ensembl_gene_id) %>%
  dplyr::arrange(desc(Weighted_Degree)) %>%
  dplyr::slice(1:1000)

# Node.Centrality.degree <- Node.Centrality %>%
#   dplyr::filter(gene %in% Human.lncRNA$ensembl_gene_id) %>%
#   dplyr::arrange(desc(degree)) %>%
#   dplyr::slice(1:1000)

# Define modules ----------------------------------------------------------

Sponge.modules <- Define_Modules(Sponge.filtered, Node.Centrality.weighted, remove.central = T, set.parallel = F)

# Module size distribution
Size.modules <- sapply(Sponge.modules, length)

# Calculate module scores ----------------------------------------------------------
BRCA.Modules.OE <- Enrichment_Modules(TCGA.expr.tumor, Sponge.modules,
                                      bin.size = 100, min.size = 10, max.size = 200, method = "OE")
METABRIC.Modules.OE <- Enrichment_Modules(METABRIC.expr, Sponge.modules,
                                          bin.size = 100, min.size = 10, max.size = 200, method = "OE")

stopCluster(cl)

BRCA.Modules.GSVA <- Enrichment_Modules(TCGA.expr.tumor, Sponge.modules, bin.size = 100, min.size = 10,
                                      max.size = 200, method = "GSVA", cores = 16)
METABRIC.Modules.GSVA <- Enrichment_Modules(METABRIC.expr, Sponge.modules, bin.size = 100, min.size = 10,
                                          max.size = 200, method = "GSVA", cores = 16)

#Find common modules

CommonModules <- intersect(rownames(BRCA.Modules.OE), rownames(METABRIC.Modules.OE))
BRCA.Modules.OE <- BRCA.Modules.OE[CommonModules, ]
METABRIC.Modules.OE <- METABRIC.Modules.OE[CommonModules, ]

# Train subtype classifier -------------------------------------------------------------
no_cores <- detectCores() - 1
cl <- makePSOCKcluster(25)
registerDoParallel(cl)

# Define input
Inputdata.model <- t(BRCA.Modules.OE) %>% scale(center = T, scale = T) %>%
  as.data.frame()
Inputdata.model <- Inputdata.model %>%
  mutate(Class = as.factor(TCGA.meta.tumor$SUBTYPE[match(rownames(Inputdata.model), TCGA.meta.tumor$sampleID)]))
Inputdata.model <- Inputdata.model[! is.na(Inputdata.model$Class), ]

# Define hyperparameters
Metric <- "Exact_match"
tunegrid <- expand.grid(.mtry=c(1:100))
n.folds <- 10 # number of folds
repetitions <- 3 # number of k-fold cv iterations

# Calibrate model
SpongingActiivty.model <- RF_Classifier(Inputdata.model, n.folds, repetitions, metric = Metric, tunegrid)

#Test classification performance on second cohort
METABRIC.Modules.OE <- METABRIC.Modules.OE[ ,complete.cases(t(METABRIC.Modules.OE))]
Meta.metabric <- METABRIC.meta$CLAUDIN_SUBTYPE[match(colnames(METABRIC.Modules.OE), METABRIC.meta$PATIENT_ID)]
Input.Metabric <- t(METABRIC.Modules.OE) %>% scale(center = T, scale = T)
Prediction.model <- predict(SpongingActiivty.model$Model, Input.Metabric)
SpongingActiivty.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.model), Meta.metabric)

stopCluster(cl)
# Define random modules and use them to train new classifier -------------------------------------------------------------
no_cores <- detectCores() - 1
cl <- makePSOCKcluster(15)
registerDoParallel(cl)

Sizes.modules <- lengths(Sponge.modules)

# Define random modules
Random.Modules <- list()

for(j in 1:length(Size.modules)) {
  Module.Elements <-  sample.int(n = nrow(TCGA.expr.tumor), size = Size.modules[j],
                                 replace = F)
  # print(Module.Elements)
  Random.Modules[[j]] <- rownames(TCGA.expr.tumor)[Module.Elements]
}
names(Random.Modules) <- names(Sponge.modules)

# Calculate enrichment scores for BRCA and METABRIC
BRCA.RandomModules.OE <- Enrichment_Modules(TCGA.expr.tumor, Random.Modules,
                                            bin.size = 100, min.size = 10, max.size = 200, method = "OE")
METABRIC.RandomModules.OE <- Enrichment_Modules(METABRIC.expr, Random.Modules,
                                                bin.size = 100, min.size = 10, max.size = 200, method = "OE")

#Find common modules

CommonModules <- intersect(rownames(BRCA.RandomModules.OE), rownames(METABRIC.RandomModules.OE))
BRCA.RandomModules.OE <- BRCA.RandomModules.OE[CommonModules, ]
METABRIC.RandomModules.OE <- METABRIC.RandomModules.OE[CommonModules, ]

# Train model
Inputdata.model.RANDOM <- t(BRCA.RandomModules.OE) %>% scale(center = T, scale = T) %>%
  as.data.frame()
Inputdata.model.RANDOM <- Inputdata.model.RANDOM %>%
  mutate(Class = as.factor(TCGA.meta.tumor$SUBTYPE[match(rownames(Inputdata.model.RANDOM), TCGA.meta.tumor$sampleID)]))

Inputdata.model.RANDOM <- Inputdata.model.RANDOM[! is.na(Inputdata.model.RANDOM$Class), ]

Random.model <- RF_Classifier(Inputdata.model.RANDOM, n.folds, repetitions, metric = Metric, tunegrid)

#Test classification performance on second cohort
METABRIC.RandomModules.OE <- METABRIC.RandomModules.OE[ ,complete.cases(t(METABRIC.RandomModules.OE))]
Meta.metabric <- METABRIC.meta$CLAUDIN_SUBTYPE[match(colnames(METABRIC.RandomModules.OE), METABRIC.meta$PATIENT_ID)]

Input.Metabric.random <- t(METABRIC.RandomModules.OE) %>% scale(center = T, scale = T)
Prediction.random.model <- predict(Random.model$Model, Input.Metabric.random)
Random.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.random.model), Meta.metabric)

stopCluster(cl)

# Build classifier only on central genes-------------------------------------------------------------
no_cores <- detectCores() - 1
cl <- makePSOCKcluster(15)
registerDoParallel(cl)

# Define central genes that are also present in test set (METABRIC)
Common.CentralGenes <- intersect(rownames(METABRIC.expr), intersect(rownames(TCGA.expr.tumor), rownames(BRCA.Modules.OE)))

# Define model input
#TCGA.expr.tumor <- TCGA.expr.tumor[, match(TCGA.meta.tumor$sampleID, colnames(TCGA.expr.tumor))]

Inputdata.centralGene <- t(TCGA.expr.tumor[rownames(TCGA.expr.tumor) %in% Common.CentralGenes, ]) %>% scale(center = TRUE, scale = TRUE) %>%
  as.data.frame()
Inputdata.centralGene <- Inputdata.centralGene %>%
  mutate(Class = as.factor(TCGA.meta.tumor$SUBTYPE[match(rownames(Inputdata.centralGene), TCGA.meta.tumor$sampleID)]))

Inputdata.centralGene <- Inputdata.centralGene[! is.na(Inputdata.centralGene$Class), ]

CentralGenes.model <- RF_Classifier(Inputdata.centralGene, n.folds, repetitions, metric = Metric, tunegrid)

#Test classification performance on second cohort
Inputdata.centralGene.Metabric <- t(METABRIC.expr[Common.CentralGenes, ]) %>% scale(center = TRUE, scale = TRUE)
Prediction.centralGenes.model <- predict(CentralGenes.model$Model, Inputdata.centralGene.Metabric)
CentralGenes.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.centralGenes.model), as.factor(METABRIC.meta$CLAUDIN_SUBTYPE))

stopCluster(cl)

save.image(file = "./Results/Objects/Workspace_lncRNAModel_13012022.RData")

