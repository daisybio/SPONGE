#' Convert hugo symbols to ensemble gene names.
#'
#' @importFrom Biobase
#' @importFrom biomaRt
#'
#' @param ceRNA_expression_data ceRNA expression data (same structure as input
#' for SPONGE)
#' (Expression matrix (genes x samples). hugo symbols are in first column)
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#' (See https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)
#' @param bioMart_gene_symbol_columns bioMart dataset column for gene symbols
#' (e.g. human: hgnc_symbol, mouse: mgi_symbol)
#'
#' @return Expression matrix (genes x samples). Row names are ensembl
#' gene symbols
Convert_Gene_Names <- function(ceRNA_expression_data,
                               bioMart_gene_ensembl,
                               bioMart_gene_symbol_columns) {

  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  not_done=TRUE

  ceRNA_expression_data <-  ceRNA_expression_data[!duplicated(ceRNA_expression_data[,1]), ]
  genes <- ceRNA_expression_data[,1] %>% as.character()

  G_list= data.frame()

  while(not_done)
  {
    tryCatch({
      ensembl <- biomaRt::useMart("ensembl", dataset=bioMart_gene_ensembl)
      GeneNames.df <- getBM(attributes=c('ensembl_gene_id',
                                         'external_gene_name'),
                            filters = bioMart_gene_symbol_columns,mart = ensembl,values = genes,
                            useCache = FALSE)
      not_done=FALSE
    }, warning = function(w) {
      print("WARNING SECTION")
      print(w)
    }, error = function(e) {
      print("ERROR SECTION")
      print(e)
    }, finally = {
    })
  }


  ceRNA_expression_data$Ensembl <- GeneNames.df[match(ceRNA_expression_data[,1], GeneNames.df$external_gene_name), 1]
  ceRNA_expression_data <- ceRNA_expression_data[!is.na(ceRNA_expression_data$Ensembl), ]

  rownames(ceRNA_expression_data) <- ceRNA_expression_data$Ensembl

  return(dplyr::select(ceRNA_expression_data[-1], -Ensembl))
}

#' Calculate z-scores
#'
#' @param x input to calculate z score on
#'
#' @return z score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#' Preprocessing ceRNA network
#' @param network  ceRNA network as data (typically present in the outputs of
#' sponge)
#' @param mscor.threshold mscor threshold (default 0.1)
#' @param padj.threshold adjusted p-value threshold (default 0.01)
#'
#' @return filtered ceRNA network
filter_network <- function(network,
                           mscor.threshold = .1,
                           padj.threshold = .01) {
  network %>%
    filter(mscor > mscor.threshold & p.adj < padj.threshold)
}

#' Function to calculate centrality scores
#' Calculation of combined centrality scores as proposed by Del Rio et al. (2009)
#'
#' @param CentralityMeasures dataframe with centrality score measures as columns
#' and samples as rows
#'
#' @return Vector containing combined centrality scores
Combined_Centrality <- function(CentralityMeasures) {

  CombinedCentrality.Score <- c()

  for (v in 1:nrow(CentralityMeasures)) {
    combined.c <- 0
    for (m in  colnames(CentralityMeasures)) {

      max.c <- max(CentralityMeasures[,m])
      min.c <- min(CentralityMeasures[,m])

      gene.c <- CentralityMeasures[v,m]

      combined.c <- combined.c + ((max.c - gene.c) / (max.c - min.c)^2)
    }

    CombinedCentrality.Score <- c(CombinedCentrality.Score, 0.5*combined.c)
  }
  return(CombinedCentrality.Score)
}

#' Function to calculate centrality scores
#' Calculation of weighted degree scores based on Opsahl et al. (2010)
#' Hyperparameter to tune: Alpha = 0 --> degree centrality as defined in
#' Freeman, 1978 (number of edges).
#' @param network Network formatted as a dataframe with three columns containing
#' respectively node1, node2 and weights
#' @param undirected directionality of the network (default: T)
#' @param Alpha degree centrality as defined in Barrat et al., 2004 (default: 1)
#'
#'
#' @return Dataframe containing information about nodes and their weighted
#' centrality measure
Weighted_degree <- function(network,
                            undirected = T,
                            Alpha = 1){

  # Format input matrix by using numeric as node IDs
  Nodes <- data.frame(Nodes = union(network$geneA, network$geneB),
                      Nodes_numeric = seq(1, length(union(network$geneA,
                                                          network$geneB))))

  geneA.numeric <- Nodes$Nodes_numeric[match(network$geneA, Nodes$Nodes)]
  geneB.numeric <- Nodes$Nodes_numeric[match(network$geneB, Nodes$Nodes)]

  Input.network <- data.frame(Sender = geneA.numeric, Receiver = geneB.numeric,
                              Weight = network$mscor)

  if (undirected) {
    # Define networks as undirected
    Undirected.net <- Input.network %>% tnet::symmetrise_w()
    Weighted_degree <- tnet::degree_w(Undirected.net, alpha = Alpha) %>%
      as.data.frame()

    Nodes$Weighted_degree <- Weighted_degree$output[match(Weighted_degree$node,
                                                          Nodes$Nodes_numeric)]
    return(Nodes)

  }
}

#' Functions to define Sponge modules, created as all the first neighbors of
#' the most central genes
#'
#' @param network Network as dataframe and list of central nodes. First two
#' columns of the dataframe should contain the information of the nodes
#' connected by edges.
#' @param  central.modules consider central gene as part of the module
#' (default: False)
#' @param remove.central Possibility of keeping or removing (default) central
#' genes in the modules (default: T)
#' @param set.parallel paralleling calculation of Define_Modules (default: T)
#'
#' @return List of modules. Module names are the corresponding central genes.
Define_Modules <- function(network,
                           central.modules = F,
                           remove.central = T,
                           set.parallel = T) {
  if (set.parallel) {
    Sponge.Modules <- foreach(i = 1:nrow(central.modules),  .packages = "dplyr") %dopar% {
      Sponge.temp <- network %>%
        dplyr::filter(geneA == central.modules$gene[i] | geneB == central.modules$gene[i])

      if(nrow(Sponge.temp) != 0) {
        Module.temp <- union(Sponge.temp$geneA, Sponge.temp$geneB)
        Module.temp <- Module.temp[Module.temp != central.modules$gene[i]]

        names(Module.temp) <- as.character(central.modules$gene[i])
        }
      }
    } else {

      Sponge.Modules <- list()
      k <- 1
      for(i in 1:nrow(central.modules)) {
        Sponge.temp <- network %>%
          filter(geneA == central.modules$gene[i] | geneB == central.modules$gene[i])

        if(nrow(Sponge.temp) != 0) {
          Module.temp <- union(Sponge.temp$geneA, Sponge.temp$geneB)
          Module.temp <- Module.temp[Module.temp != central.modules$gene[i]]

          Sponge.Modules[[k]] <- Module.temp
          names(Sponge.Modules)[k] <- as.character(central.modules$gene)[i]

          k <- k+1

        }
      }
    }

  if (remove.central) {
    # Remove central lncRNA from each module
    Sponge.Modules.DoubleRemoved <- list()

    for (z in 1:length(Sponge.Modules)) {
      toKeep <- !(Sponge.Modules[[z]] %in% central.modules$gene)
      Sponge.Modules.DoubleRemoved[[z]] <- Sponge.Modules[[z]][!(Sponge.Modules[[z]] %in% central.modules$gene)]
    }

    names(Sponge.Modules.DoubleRemoved) <- names(Sponge.Modules)
    return(Sponge.Modules.DoubleRemoved)
  } else {return(Sponge.Modules)}
}

#' Function to calculate enrichment scores of modules OE
#' (functions taken from: Jerby-Arnon et al. 2018)
#' @param NormCount normalized counts
#' @param gene.sign significant genes
#' @param bin.size bin size (default: 100)
#' @param num.rounds number of rounds (default: 1000)
#' @param set_seed seed size (default: 42)
#'
#' @return Signature scores
OE_module <- function(NormCount,
                      gene.sign,
                      bin.size = 100,
                      num.rounds = 1000,
                      set_seed = 42){
  set.seed(set_seed)

  SequencedGenes <- rownames(NormCount)
  # Center gene expression
  Centered.expr <- NormCount %>% t() %>% scale(center =T, scale = F) %>% t()

  # Average expression of a gene across cells
  Genes.dist <- rowMeans(NormCount,na.rm = T)

  # Divide genes in  expression bins
  Discretized.Genes.dist <- discretize(Genes.dist,n.cat = bin.size)

  # Calculate OE scores
  Signature.scores <- matrix(data = 0,nrow = ncol(NormCount),ncol =1)
  sign.names <- names(gene.sign)
  colnames(Signature.scores) <- sign.names

  Signature.scores.raw <- Signature.scores
  random.scores<- Signature.scores

  b.sign <- is.element(SequencedGenes, gene.sign)
  if(sum(b.sign)<2) {next()}

  rand.scores<-get.semi.random.OE(Centered.expr, Discretized.Genes.dist,b.sign,num.rounds = num.rounds)

  raw.scores <- colMeans(Centered.expr[b.sign,])
  final.scores <- raw.scores-rand.scores
  Signature.scores[,1] <- final.scores
  Signature.scores.raw[,1] <- raw.scores
  random.scores[,1] <- rand.scores

  return(Signature.scores)
}

#' Function to calculate semi random enrichment scores of modules OE
#' (functions taken from: Jerby-Arnon et al. 2018)
#' @param r expression matrix
#' @param genes.dist.q values of the genes after binning (result of binning)
#' @param b.sign does the signature contain less than 2 genes?
#' (controll parameter) (is set by mother function (OE module function))
#' @param num.rounds number of rounds (default: 1000)
#'
#' @return random signature scores
get.semi.random.OE <- function(r,
                               genes.dist.q,
                               b.sign,
                               num.rounds = 1000){
  # Previous name: get.random.sig.scores
  sign.q<-as.matrix(table(genes.dist.q[b.sign]))
  q<-rownames(sign.q)
  idx.all<-c()
  B<-matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
  Q<-matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)
  for (i in 1:nrow(sign.q)){
    num.genes<-sign.q[i]
    if(num.genes>0){
      idx<-which(is.element(genes.dist.q,q[i]))
      for (j in 1:num.rounds){
        idxj<-sample(idx,num.genes)
        Q[i,j]<-sum(B[idxj,j]==T)
        B[idxj,j]<-T
      }
    }
  }
  rand.scores<-apply(B,2,function(x) colMeans(r[x,]))

  rand.scores<-rowMeans(rand.scores)
  return(rand.scores)
}

#' discretize
#' @param v gene distance (defined by mother function OE module function)
#' @param n.cat size of the bins (defined by mother function  OE module function)
#'
#' @return discreted
discretize<-function(v,
                     n.cat){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
  u<-matrix(nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<q1[i])]<-i
  }
  return(u)
}

# 3) Calculate enrichment scores
# Input: List of modules (with module names), expression data (gene x samples)
# Output: matrix containing module enrichment scores (module x samples)

#' Calculate enrichment scores
#'
#' @importFrom tidyverse
#' @importFrom caret
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom dplyr
#'
#' @param Expr.matrix ceRNA expression matrix
#' @param modules Result of Define_Modules
#' @param bin.size bin size (default: 100)
#' @param min.size minimum module size (default: 10)
#' @param max.size maximum module size (default: 200)
#' @param method Enrichment to be used (Overall Enrichment: OE or Gene Set
#' Variation Analysis: GSVA) (default: OE)
#' @param cores cores for parallelisation (default: 1)
#'
#' @return matrix containing module enrichment scores (module x samples)
Enrichment_Modules <- function(Expr.matrix,
                               modules,
                               bin.size = 100,
                               min.size = 10,
                               max.size = 200,
                               method = "OE",
                               cores = 1) {

  print(paste0("Calculating modules with bin size: ", bin.size, ", min size: ", min.size, ", max size:", max.size))

  if (method == "OE") {
    Enrichmentscores.modules <- foreach(Module = 1:length(modules), .packages = c("dplyr", "GSVA"),
                                        .export = c("OE_module", "get.semi.random.OE", "discretize"),
                                        .combine = "cbind") %dopar% {
      results <- list()

      if(length(modules[[Module]]) > min.size & length(modules[[Module]]) < max.size) {
        if (sum((modules[[Module]] %in% rownames(Expr.matrix)), na.rm = TRUE) > 10) {
          Enrichment.module <- round(OE_module(Expr.matrix,modules[[Module]],bin.size),2)
          colnames(Enrichment.module) <- names(modules)[Module]
          return(Enrichment.module)
        }
        }
      }
    Enrichmentscores.modules <- Enrichmentscores.modules %>% t() %>% as.data.frame

    colnames(Enrichmentscores.modules) <- colnames(Expr.matrix)

  } else if (method == "GSVA") {
    Enrichmentscores.modules <- GSVA::gsva(as.matrix(Expr.matrix), modules, min.sz= min.size,
                                     max.sz=200,method= 'gsva',parallel.sz = cores,
                                     verbose=FALSE) %>% as.data.frame()

    colnames(Enrichmentscores.modules) <- colnames(Expr.matrix)

  }

  if (!is_empty(Enrichmentscores.modules)) {
  	return(Enrichmentscores.modules)
  } else {
    return(NULL)
  }
}

#' Calibrate classification method
#' @param data Dataframe with module scores/covariates (modules x samples)
#' AND outcome variable
#' @param lev (default: NULL)
#' @param model (default: NULL)
#'
#' @return Model and confusion matrix in a list
ExactMatch_Summary <- function (data,
                                lev = NULL,
                                model = NULL) {
  metric <- sum(data$obs == data$pred, na.rm = TRUE) / length(data$obs)
  names(metric) <- "Exact_match"
  return(metric)
}

#' RF classification model
#'
#' @param Input.object data.frame made by predictors and dependent variable
#' @param K number of folds (k-fold)
#' @param rep number of times repeating the cross validation
#' @param metric metric (Exact_match, Accuracy) (default: Exact_match)
#' @param tunegrid defines the grid for the hyperparameter optimization during
#' cross validation (caret package)
#' @param set_seed set seed (default: 42)
#'
#' @return
RF_Classifier <- function(Input.object,
                          K,
                          rep,
                          metric = "Exact_match",
                          tunegrid,
                          set_seed = 42){
  set.seed(set_seed)

  # Training settings
  # k Fold stratified CV
  trainIndex <- createFolds(Input.object$Class, k = K,
                            list = T, returnTrain = T)

  if (metric == "Exact_match") {
    control <- trainControl(method="repeatedcv", number=K, repeats=rep, search="grid",
                            savePredictions = "final", classProbs = F,
                            summaryFunction = ExactMatch_Summary,
                            allowParallel = TRUE, index = trainIndex)
  } else if (metric == "Accuracy") {
      control <- trainControl(method="repeatedcv", number=K, repeats=rep, search="grid",
                              savePredictions = "final", classProbs = F,
                              allowParallel = TRUE, index = trainIndex)
  }

  rf.model <- train(Class~., data=Input.object, method="rf", metric=metric,
                    tuneGrid=tunegrid, trControl=control,
                    importance = TRUE)

  Confusion.matrix <- confusionMatrix(rf.model$pred$pred, rf.model$pred$obs)

  Output.rf <-list(Model = rf.model, ConfusionMatrix_training = Confusion.matrix)
  return(Output.rf)
}

#' prepare TCGA formats for spongEffects
#'
#' @importFrom tidyverse
#' @importFrom caret
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom dplyr
#'
#' @param tcga_cancer_symbol e.g., BRCA for breast cancer
#' @param normal_ceRNA_expression_data normal ceRNA expression data
#' (same structure as input for SPONGE)
#' @param tumor_ceRNA_expression_data tumor ceRNA expression data
#' (same structure as input for SPONGE)
#' @param normal_metadata metadata for normal samples
#' (TCGA format style, needs to include column: sampleID, PATIENT_ID)
#' @param tumor_metadata metadata for tumor samples
#' (TCGA format style, needs to include column: sampleID, PATIENT_ID)
#' @param clinical_data clinical data for all patients
#' (TCGA format style, needs to include column: PATIENT_ID, AJCC_PATHOLOGIC_TUMOR_STAGE)
#' @param tumor_stages_of_interest array e.g.,
#' c(STAGE I', 'STAGE IA', 'STAGE IB', 'STAGE II', 'STAGE IIA')
#' @param subtypes_of_interest array e.g.,
#' c("LumA",  "LumB",  "Her2",  "Basal", "Normal")
#'
#' @return list of prepared data. You can access it with list$objectname for
#' further spongEffects steps
prepare_tcga_for_spongEffects <- function(tcga_cancer_symbol,
                                          normal_ceRNA_expression_data,
                                          tumor_ceRNA_expression_data,
                                          normal_metadata,
                                          tumor_metadata,
                                          clinical_data,
                                          tumor_stages_of_interest,
                                          subtypes_of_interest){

  tcga_cancer_symbol=paste0(tcga_cancer_symbol,"_")

  # Tidy and prepare input dataset
  TCGA.expr.tumor <- tumor_ceRNA_expression_data %>% t() %>% as.data.frame()
  TCGA.expr.normal <-normal_ceRNA_expression_data %>% t() %>% as.data.frame()
  TCGA.meta.normal <- normal_metadata %>%
    dplyr::filter(sampleID %in% colnames(TCGA.expr.normal))
  TCGA.meta.tumor <- tumor_metadata %>%
    dplyr::filter(sampleID %in% colnames(TCGA.expr.tumor))

  TCGA.meta.tumor$PATIENT_ID <- substr(TCGA.meta.tumor$sampleID,1,nchar(TCGA.meta.tumor$sampleID)-3)

  # Load metadata with grading information and keep only samples with StageI, II, III, and IV
  # clean stage mess
  MetaData_TCGA_Complete <- clinical_data %>%
    dplyr::filter(!(SUBTYPE %in% c("")) & PATIENT_ID %in% TCGA.meta.tumor$PATIENT_ID) %>%
    dplyr::filter(AJCC_PATHOLOGIC_TUMOR_STAGE %in% tumor_stages_of_interest) %>%
    mutate(Grading = ifelse(AJCC_PATHOLOGIC_TUMOR_STAGE %in% c('STAGE I', 'STAGE IA', 'STAGE IB'), "Stage_I",
                            ifelse(AJCC_PATHOLOGIC_TUMOR_STAGE %in% c('STAGE II', 'STAGE IIA', 'STAGE IIB'), "Stage_II",
                                   ifelse(AJCC_PATHOLOGIC_TUMOR_STAGE %in% "STAGE IV", "Stage_IV", "Stage_III"))))

  TCGA.meta.tumor <- TCGA.meta.tumor[match(MetaData_TCGA_Complete$PATIENT_ID, TCGA.meta.tumor$PATIENT_ID), ] %>%
    left_join(MetaData_TCGA_Complete, .by = PATIENT_ID)
  TCGA.meta.tumor$SUBTYPE <- gsub(tcga_cancer_symbol,"", TCGA.meta.tumor$SUBTYPE)
  TCGA.meta.tumor$SUBTYPE <- factor(TCGA.meta.tumor$SUBTYPE, levels = subtypes_of_interest)

  TCGA.expr.normal <- TCGA.expr.normal[, match(TCGA.meta.normal$sampleID, colnames(TCGA.expr.normal))]
  TCGA.expr.tumor <- TCGA.expr.tumor[, match(TCGA.meta.tumor$sampleID, colnames(TCGA.expr.tumor))]

  rownames(TCGA.expr.tumor) <- gsub("\\..*","",rownames(TCGA.expr.tumor))
  rownames(TCGA.expr.normal) <- gsub("\\..*","",rownames(TCGA.expr.normal))

  prep_TCGA_spongEffects <- list(TCGA.expr.tumor,TCGA.expr.normal,TCGA.meta.tumor,TCGA.meta.normal,MetaData_TCGA_Complete)
  names(prep_TCGA_spongEffects) <- c("TCGA.expr.tumor","TCGA.expr.normal","TCGA.meta.tumor","TCGA.meta.normal","MetaData_TCGA_Complete")

  return(prep_TCGA_spongEffects)
}
#' prepare METABRIC formats for spongEffects
#'
#' @importFrom tidyverse
#' @importFrom caret
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom dplyr
#' @importFrom Biobase
#'
#' @param metabric_expression filepath to expression data in metabric format
#' @param metabric_metadata filepath to metabric metadata in metabric format
#' @param subtypes_of_interest array e.g.,
#' c("LumA",  "LumB",  "Her2",  "Basal", "Normal")
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#' (See https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)
#' (default: hsapiens_gene_ensembl)
#' @param bioMart_gene_symbol_columns bioMart dataset column for gene symbols
#' (e.g. human: hgnc_symbol, mouse: mgi_symbol)
#' (default: hgnc_symbol)
#'
#' @return list with metabric expression and metadata. You can access it with
#' list$objectname for further spongEffects steps
prepare_metabric_for_spongEffects <- function(metabric_expression,
                                              metabric_metadata,
                                              subtypes_of_interest,
                                              bioMart_gene_ensembl = "hsapiens_gene_ensembl",
                                              bioMart_gene_symbol_columns = "hgnc_symbol"){
  METABRIC.expr <- read.delim(metabric_expression,
                              header=T) %>% dplyr::select(-Entrez_Gene_Id)
  METABRIC.expr <- METABRIC.expr[Biobase::isUnique(METABRIC.expr$Hugo_Symbol), ]

  METABRIC.expr <- Convert_Gene_Names(METABRIC.expr,bioMart_gene_ensembl,bioMart_gene_symbol_columns)
  colnames(METABRIC.expr) <- gsub("\\.", "-", colnames(METABRIC.expr))

  METABRIC.meta <- read.delim(metabric_metadata, comment.char="#") %>%
    filter(CLAUDIN_SUBTYPE %in% subtypes_of_interest & PATIENT_ID %in% colnames(METABRIC.expr))
  METABRIC.meta$CLAUDIN_SUBTYPE <- factor(METABRIC.meta$CLAUDIN_SUBTYPE, levels = subtypes_of_interest)

  METABRIC.expr <- METABRIC.expr[, match(METABRIC.meta$PATIENT_ID, colnames(METABRIC.expr))]

  prep_metabric_spongEffects<-list(METABRIC.expr,METABRIC.meta)
  names(prep_metabric_spongEffects) <- c("METABRIC.expr","METABRIC.meta")

  return(prep_metabric_spongEffects)
}

#' prepare ceRNA network and network centralities from SPONGE / SPONGEdb
#' for spongEffects
#'
#' @importFrom tidyverse
#' @importFrom caret
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom dplyr
#'
#' @param sponge_effects the ceRNA network downloaded as R object from SPONGEdb
#' (Hoffmann et al., 2021) or created by SPONGE (List et al., 2019)
#' usually ends with _sponge_results
#' @param network_analysis the network analysis downloaded as R object
#' from SPONGEdb (Hoffmann et al., 2021) or created by SPONGE
#' (List et al., 2019) usually ends with _networkAnalysis
#' if network_analysis is NA then the function only filters the ceRNA network
#' @param mscor.threshold mscor threshold to be filtered (default: 0.1)
#' @param padj.threshold adjusted p-value to be filtered (default: 0.01)
#'
#' @return list of filtered ceRNA network and network centrailies. You can
#' access it with list$objectname for further spongEffects steps
filter_ceRNA_network <- function(sponge_effects,
                                 network_analysis = NA,
                                 mscor.threshold = 0.1,
                                 padj.threshold = 0.01){

  #Filter SPONGE network for significant edges
  Sponge.filtered <- sponge_effects %>%
    filter_network(mscor.threshold =  .1, padj.threshold = .01)

  # Calculate weighted centrality scores and add them to the ones present in SpongeDB
  #cancerType.Centrality <- gsub(" ", "_", cancerType)

  if(is.na(network_analysis))
  {
    sponge_network_centralites <- list(Sponge.filtered)
    names(sponge_network_centralites) <- c("Sponge.filtered")
  }
  else
  {
    Node.Centrality <- network_analysis
    Nodes <- Weighted_degree(Sponge.filtered, undirected = T, Alpha = 1)

    Node.Centrality <- Node.Centrality %>%
      mutate(Weighted_Degree = Nodes$Weighted_degree[match(Node.Centrality$gene, Nodes$Nodes)])
    sponge_network_centralites <- list(Sponge.filtered,Node.Centrality)
    names(sponge_network_centralites) <- c("Sponge.filtered","Node.Centrality")

  }
  return(sponge_network_centralites)

}


#' prepare ceRNA network and network centralities from SPONGE / SPONGEdb
#' @importFrom tidyverse
#' @importFrom caret
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom dplyr
#' @importFrom Biobase
#' @importFrom biomaRt
#'
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#' (See https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)
#' @param weighted_node_centrality output from filter_ceRNA_network()
#' @param cutoff the top cutoff modules will be returned (default: 1000)
#'
#' @return top cutoff modules, with lncRNA as central genes
get_lncRNA_modules <- function(bioMart_gene_ensembl,
                               weighted_node_centrality,
                               cutoff = 1000){
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  not_done=TRUE
  while(not_done)
  {
    tryCatch({
      ensembl=useMart("ensembl")
      ensemblH = useDataset(bioMart_gene_ensembl,mart=ensembl)
      Human.lncRNA = getBM(attributes=c("ensembl_gene_id","gene_biotype","description"), mart=ensemblH, useCache = FALSE) %>%
        dplyr::filter(gene_biotype == "lncRNA")

      not_done=FALSE
    }, warning = function(w) {
      print("WARNING SECTION")
      print(w)
    }, error = function(e) {
      print("ERROR SECTION")
      print(e)
    }, finally = {
    })
  }


  Node.Centrality.weighted <- weighted_node_centrality %>%
    dplyr::filter(gene %in% Human.lncRNA$ensembl_gene_id) %>%
    dplyr::arrange(desc(Weighted_Degree)) %>%
    dplyr::slice(1:cutoff)

  return(Node.Centrality.weighted)
}






