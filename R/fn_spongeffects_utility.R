#' Convert hugo symbols to ensemble gene names.
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
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
fn_convert_gene_names <- function(ceRNA_expression_data,
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

#' #' Calculate z-scores
#' #'
#' #' @import tidyverse
#' #' @import caret
#' #' @import dplyr
#' #' @import Biobase
#' #' @import biomaRt
#' #' @import randomForest
#' #' @import ggridges
#' #' @import cvms
#' #' @import miRBaseConverter
#' #'
#' #' @param x input to calculate z score on
#' #'
#' #' @return z score
#' fn_cal_z_score <- function(x){
#'   (x - mean(x)) / sd(x)
#' }

#' Preprocessing ceRNA network
#' @param network  ceRNA network as data (typically present in the outputs of
#' sponge)
#' @param mscor.threshold mscor threshold (default 0.1)
#' @param padj.threshold adjusted p-value threshold (default 0.01)
#'
#' @return filtered ceRNA network
fn_filter_network <- function(network,
                              mscor.threshold = .1,
                              padj.threshold = .01) {
    network %>%
        filter(mscor > mscor.threshold & p.adj < padj.threshold)
}

#' Function to calculate centrality scores
#' Calculation of combined centrality scores as proposed by Del Rio et al. (2009)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param CentralityMeasures dataframe with centrality score measures as columns
#' and samples as rows
#'
#' @return Vector containing combined centrality scores
fn_combined_centrality <- function(CentralityMeasures) {

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
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import tnet
#'
#' @param network Network formatted as a dataframe with three columns containing
#' respectively node1, node2 and weights
#' @param undirected directionality of the network (default: T)
#' @param Alpha degree centrality as defined in Barrat et al., 2004 (default: 1)
#'
#'
#' @return Dataframe containing information about nodes and their weighted
#' centrality measure
fn_weighted_degree <- function(network,
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
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param network Network as dataframe and list of central nodes. First two
#' columns of the dataframe should contain the information of the nodes
#' connected by edges.
#' @param  central.modules consider central gene as part of the module
#' (default: False)
#' @param remove.central Possibility of keeping or removing (default) central
#' genes in the modules (default: T)
#' @param set.parallel paralleling calculation of define_modules() (default: F)
#'
#' @export
#'
#' @return List of modules. Module names are the corresponding central genes.
define_modules <- function(network,
                           central.modules = F,
                           remove.central = T,
                           set.parallel = T) {
    set.parallel = F
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

#' discretize
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param v gene distance (defined by mother function OE module function)
#' @param n.cat size of the bins (defined by mother function  OE module function)
#'
#' @return discretized
fn_discretize_spongeffects<-function(v,
                                     n.cat){
    q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
    u<-matrix(nrow = length(v))
    for(i in 2:n.cat){
        u[(v>=q1[i-1])&(v<q1[i])]<-i
    }
    return(u)
}

#' Function to calculate enrichment scores of modules OE
#' (functions taken from: Jerby-Arnon et al. 2018)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param NormCount normalized counts
#' @param gene.sign significant genes
#' @param bin.size bin size (default: 100)
#' @param num.rounds number of rounds (default: 1000)
#' @param set_seed seed size (default: 42)
#'
#' @return Signature scores
fn_OE_module <- function(NormCount,
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
    Discretized.Genes.dist <- fn_discretize_spongeffects(Genes.dist,n.cat = bin.size)

    # Calculate OE scores
    Signature.scores <- matrix(data = 0,nrow = ncol(NormCount),ncol =1)
    sign.names <- names(gene.sign)
    colnames(Signature.scores) <- sign.names

    Signature.scores.raw <- Signature.scores
    random.scores<- Signature.scores

    b.sign <- is.element(SequencedGenes, gene.sign)
    if(sum(b.sign)<2) {next()}

    rand.scores<-fn_get_semi_random_OE(Centered.expr, Discretized.Genes.dist,b.sign,num.rounds = num.rounds)

    raw.scores <- colMeans(Centered.expr[b.sign,])
    final.scores <- raw.scores-rand.scores
    Signature.scores[,1] <- final.scores
    Signature.scores.raw[,1] <- raw.scores
    random.scores[,1] <- rand.scores

    return(Signature.scores)
}

#' Function to calculate semi random enrichment scores of modules OE
#' (functions taken from: Jerby-Arnon et al. 2018)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param r expression matrix
#' @param genes.dist.q values of the genes after binning (result of binning)
#' @param b.sign does the signature contain less than 2 genes?
#' (controll parameter) (is set by mother function (OE module function))
#' @param num.rounds number of rounds (default: 1000)
#'
#' @return random signature scores
fn_get_semi_random_OE <- function(r,
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

    rand.scores<-apply(B,2,function(x){
        sel_matrix <- r[x,]
        if(is.array(sel_matrix)) return(colMeans(sel_matrix))
        else return(mean(sel_matrix))
    })
    if(is.array(rand.scores)) rand.scores<-rowMeans(rand.scores)
    else rand.scores <- mean(rand.scores)

    return(rand.scores)
}


# 3) Calculate enrichment scores
# Input: List of modules (with module names), expression data (gene x samples)
# Output: matrix containing module enrichment scores (module x samples)

#' Calculate enrichment scores
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import rlang
#'
#' @param Expr.matrix ceRNA expression matrix
#' @param modules Result of define_modules()
#' @param bin.size bin size (default: 100)
#' @param min.size minimum module size (default: 10)
#' @param max.size maximum module size (default: 200)
#' @param min.expr minimum expression (default: 10)
#' @param method Enrichment to be used (Overall Enrichment: OE or Gene Set
#' Variation Analysis: GSVA) (default: OE)
#' @param cores number of cores to be used to calculate entichment scores
#' with gsva or ssgsea methods. Default 1
#' @export
#'
#' @return matrix containing module enrichment scores (module x samples)
enrichment_modules <- function(gene_expr,
                               modules,
                               bin.size = 100,
                               min.size = 10,
                               max.size = 200,
                               min.expr = 10,
                               method = "OE",
                               cores = 1) {

    Expr.matrix <- gene_expr

    print(paste0("Calculating modules with bin size: ", bin.size, ", min size: ", min.size, ", max size:", max.size))

    if (method == "OE") {
        Enrichmentscores.modules <- foreach(Module = 1:length(modules), .packages = c("dplyr", "GSVA"),
                                            .export = c("fn_OE_module", "fn_get_semi_random_OE", "fn_discretize_spongeffects"),
                                            .combine = "cbind") %dopar% {
                                                results <- list()

                                                if(length(modules[[Module]]) > min.size & length(modules[[Module]]) < max.size) {
                                                    if (sum((modules[[Module]] %in% rownames(Expr.matrix)), na.rm = TRUE) > min.expr) {
                                                        Enrichment.module <- round(fn_OE_module(Expr.matrix,modules[[Module]],bin.size),2)
                                                        colnames(Enrichment.module) <- names(modules)[Module]
                                                        return(Enrichment.module)
                                                    }
                                                }
                                            }
        Enrichmentscores.modules <- Enrichmentscores.modules %>% t() %>% as.data.frame

        colnames(Enrichmentscores.modules) <- colnames(Expr.matrix)

    } else {
        Enrichmentscores.modules <- GSVA::gsva(as.matrix(Expr.matrix), modules, min.sz= min.size,
                                               max.sz=max.size,method= method,parallel.sz = cores,
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
fn_exact_match_summary <- function (data,
                                    lev = NULL,
                                    model = NULL) {
    metric <- sum(data$obs == data$pred, na.rm = TRUE) / length(data$obs)
    names(metric) <- "Exact_match"
    return(metric)
}

#' RF classification model
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
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
fn_RF_classifier <- function(Input.object,
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
                                summaryFunction = fn_exact_match_summary,
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
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
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
#' @export
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
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
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
#' @export
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

    METABRIC.expr <- fn_convert_gene_names(METABRIC.expr,bioMart_gene_ensembl,bioMart_gene_symbol_columns)
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
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param sponge_effects the ceRNA network downloaded as R object from SPONGEdb
#' (Hoffmann et al., 2021) or created by SPONGE (List et al., 2019)
#' usually ends with _sponge_results
#' @param network_analysis the network analysis downloaded as R object
#' from SPONGEdb (Hoffmann et al., 2021) or created by SPONGE
#' (List et al., 2019) usually ends with _networkAnalysis
#' if network_analysis is NA then the function only filters the ceRNA network
#' @param mscor.threshold mscor threshold to be filtered (default: NA)
#' @param padj.threshold adjusted p-value to be filtered (default: NA)
#'
#' @export
#'
#' @return list of filtered ceRNA network and network centrailies. You can
#' access it with list$objectname for further spongEffects steps
filter_ceRNA_network <- function(sponge_effects,
                                 network_analysis = NA,
                                 mscor.threshold = 0.1,
                                 padj.threshold = 0.01){

    #Filter SPONGE network for significant edges
    Sponge.filtered <- sponge_effects %>%
        fn_filter_network(mscor.threshold =  mscor.threshold, padj.threshold = padj.threshold)



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
        Nodes <- fn_weighted_degree(Sponge.filtered, undirected = T, Alpha = 1)

        Node.Centrality <- Node.Centrality %>%
            mutate(Weighted_Degree = Nodes$Weighted_degree[match(Node.Centrality$gene, Nodes$Nodes)])
        sponge_network_centralites <- list(Sponge.filtered,Node.Centrality)
        names(sponge_network_centralites) <- c("Sponge.filtered","Node.Centrality")

    }
    return(sponge_network_centralites)

}

#' prepare ceRNA network and network centralities from SPONGE / SPONGEdb
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import tnet
#'
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#' (See https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)
#' @param weighted_node_centrality output from filter_ceRNA_network()
#' @param ceRNA_class default c("lncRNA","circRNA","protein_coding") (see http://www.ensembl.org/info/genome/genebuild/biotypes.html)
#' @param cutoff the top cutoff modules will be returned (default: 1000)
#'
#' @export
#'
#' @return top cutoff modules, with lncRNA as central genes
get_central_modules <- function(bioMart_gene_ensembl,
                                weighted_node_centrality,
                                ceRNA_class = c("lncRNA","circRNA","protein_coding"),
                                cutoff = 1000){
    httr::set_config(httr::config(ssl_verifypeer = FALSE))
    not_done=TRUE
    while(not_done)
    {
        tryCatch({
            ensembl=useMart("ensembl")
            ensemblH = useDataset(bioMart_gene_ensembl,mart=ensembl)
            Human.lncRNA = getBM(attributes=c("ensembl_gene_id","gene_biotype","description"), mart=ensemblH, useCache = FALSE)

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

    all_ceRNAs <- Human.lncRNA
    Human.lncRNA<-Human.lncRNA%>% dplyr::filter(gene_biotype %in% ceRNA_class)

    c_interesting_ceRNAs<-Human.lncRNA$ensembl_gene_id

    if("circRNA" %in% ceRNA_class)
    {
        df_circRNAS_to_add <- weighted_node_centrality[ with(weighted_node_centrality, grepl("circ", gene) | grepl(":", gene) | grepl("chr", gene))]
        c_interesting_ceRNAs<- c(c_interesting_ceRNAs,df_circRNAS_to_add$gene)
    }

    Node.Centrality.weighted <- weighted_node_centrality %>%
        dplyr::filter(gene %in% c_interesting_ceRNAs) %>%
        dplyr::arrange(desc(Weighted_Degree)) %>%
        dplyr::slice(1:cutoff)

    return(Node.Centrality.weighted)
}

#' tests and trains a model for a disease using a training and test data set
#' (e.g., TCGA-BRCA and METABRIC)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param modules return from enrichment_modules() function
#' @param modules_metadata metadata table containing information about samples/patients
#' @param label Column of metadata to use as label in classification model
#' @param sampleIDs Column of metadata containing sample/patient IDs to be matched with column names of spongEffects scores
#' @param Metric metric (Exact_match, Accuracy) (default: Exact_match)
#' @param tunegrid_c defines the grid for the hyperparameter optimization during
#' cross validation (caret package) (default: 1:100)
#' @param n_folds number of folds (default: 10)
#' @param repetitions  number of k-fold cv iterations (default: 3)
#'
#' @export
#'
#' @return returns a list with the trained model and the prediction results
# train_and_test_model <- function(train_modules,
#                                  train_modules_metadata,
#                                  test_modules,
#                                  test_modules_metadata,
#                                  test_modules_meta_data_type = "TCGA",
#                                  metric = "Exact_match",
#                                  tunegrid_c = c(1:100),
#                                  n_folds = 10,
#                                  repetitions = 3){
#
#     BRCA.Modules.OE <- train_modules
#     METABRIC.Modules.OE <-test_modules
#     TCGA.meta.tumor = train_modules_metadata
#     METABRIC.meta = test_modules_metadata
#     Metric <- metric
#
#     #Find common modules
#
#     CommonModules <- intersect(rownames(BRCA.Modules.OE), rownames(METABRIC.Modules.OE))
#     BRCA.Modules.OE <- BRCA.Modules.OE[CommonModules, ]
#     METABRIC.Modules.OE <- METABRIC.Modules.OE[CommonModules, ]
#
#     # Train subtype classifier -------------------------------------------------------------
#
#     # Define input
#     Inputdata.model <- t(BRCA.Modules.OE) %>% scale(center = T, scale = T) %>%
#         as.data.frame()
#     Inputdata.model <- Inputdata.model %>%
#         mutate(Class = as.factor(TCGA.meta.tumor$SUBTYPE[match(rownames(Inputdata.model), TCGA.meta.tumor$sampleID)]))
#     Inputdata.model <- Inputdata.model[! is.na(Inputdata.model$Class), ]
#
#     # Define hyperparameters
#     Metric <- Metric
#     tunegrid <- expand.grid(.mtry=tunegrid_c)
#     n.folds <- n_folds # number of folds
#     repetitions <- repetitions # number of k-fold cv iterations
#
#     # Calibrate model
#     SpongingActivity.model <- fn_RF_classifier(Inputdata.model, n.folds, repetitions, metric = Metric, tunegrid)
#
#     #Test classification performance on second cohort
#     METABRIC.Modules.OE <- METABRIC.Modules.OE[ ,complete.cases(t(METABRIC.Modules.OE))]
#     if(test_modules_meta_data_type == "METABRIC"){
#         Meta.metabric <- METABRIC.meta$CLAUDIN_SUBTYPE[match(colnames(METABRIC.Modules.OE), METABRIC.meta$PATIENT_ID)]
#     }
#     if(test_modules_meta_data_type == "TCGA"){
#         Meta.metabric <- METABRIC.meta$SUBTYPE[match(colnames(METABRIC.Modules.OE), METABRIC.meta$sampleID)]
#     }
#
#     Input.Metabric <- t(METABRIC.Modules.OE) %>% scale(center = T, scale = T)
#     Prediction.model <- predict(SpongingActivity.model$Model, Input.Metabric)
#     Meta.metabric<-as.factor(Meta.metabric)
#     SpongingActivity.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.model), Meta.metabric)
#
#     prediction_model<-list(SpongingActivity.model,Prediction.model)
#     names(prediction_model) <- c("SpongingActivity.model","Prediction.model")
#
#     return(prediction_model)
# }

#' Calibrate classification RF classification model
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param modules return from enrichment_modules() function
#' @param modules_metadata metadata table containing information about samples/patients
#' @param label Column of metadata to use as label in classification model
#' @param sampleIDs Column of metadata containing sample/patient IDs to be matched with column names of spongEffects scores
#' @param Metric metric (Exact_match, Accuracy) (default: Exact_match)
#' @param tunegrid_c defines the grid for the hyperparameter optimization during
#' cross validation (caret package) (default: 1:100)
#' @param n_folds number of folds (default: 10)
#' @param repetitions  number of k-fold cv iterations (default: 3)
#'
#' @export
#'
#' @return returns a list with the trained model and the prediction results
calibrate_model <- function(modules,
                            modules_metadata,
                            label,
                            sampleIDs,
                            Metric = "Exact_match",
                            tunegrid_c = c(1:100),
                            n_folds = 10,
                            repetitions = 3){

    if (label %in% colnames(modules_metadata) & sampleIDs %in% colnames(modules_metadata)) {

        # Train subtype classifier -------------------------------------------------------------
        # Define input
        Inputdata.model <- t(modules) %>% scale(center = T, scale = T) %>%
            as.data.frame()
        Inputdata.model <- Inputdata.model %>%
            mutate(Class = as.factor(modules_metadata[match(rownames(Inputdata.model), TCGA.meta.tumor[,sampleIDs]), label]))
        Inputdata.model <- Inputdata.model[! is.na(Inputdata.model$Class), ]

        # Define hyperparameters
        Metric <- Metric
        tunegrid <- expand.grid(.mtry=tunegrid_c)
        n.folds <- n_folds # number of folds
        repetitions <- repetitions # number of k-fold cv iterations

        # Calibrate model
        SpongingActivity.model <- fn_RF_classifier(Inputdata.model, n.folds, repetitions, metric = Metric, tunegrid)

        return(prediction_model)

    } else {print("label and/or sampleIDs must be columns in metadata")}
}

#' build classifiers for central genes
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param train_gene_expr expression data of train dataset,
#' genenames must be in rownames
#' @param test_gene_expr expression data of test dataset,
#' genenames must be in rownames
#' @param train_enrichment_modules return of enrichment_modules()
#' @param test_enrichment_modules return of enrichment_modules()
#' @param train_meta_data meta data of train dataset
#' @param test_meta_data meta data of test dataset
#' @param train_meta_data_type TCGA or METABRIC
#' @param test_meta_data_type TCGA or METABRIC
#' @param metric metric (Exact_match, Accuracy) (default: Exact_match)
#' @param tunegrid_c defines the grid for the hyperparameter optimization during
#' cross validation (caret package) (default: 1:100)
#' @param n.folds number of folds to be calculated
#' @param repetitions  number of k-fold cv iterations (default: 3)
#'
#' @export
#'
#' @return model for central genes
build_classifier_central_genes<-function(train_gene_expr,
                                         test_gene_expr,
                                         train_enrichment_modules,
                                         test_enrichment_modules,
                                         train_meta_data,
                                         test_meta_data,
                                         train_meta_data_type="TCGA",
                                         test_meta_data_type="TCGA",
                                         metric="Exact_match",
                                         tunegrid_c=c(1:100),
                                         n.folds = 10,
                                         repetitions=3){

    TCGA.expr.tumor<-train_gene_expr
    METABRIC.expr<-test_gene_expr
    BRCA.Modules.OE<-train_enrichment_modules
    TCGA.meta.tumor<-train_meta_data
    METABRIC.meta<-test_meta_data
    Metric<-metric

    tunegrid<-expand.grid(.mtry=tunegrid_c)


    # Define central genes that are also present in test set (METABRIC)
    Common.CentralGenes <- intersect(rownames(METABRIC.expr), intersect(rownames(TCGA.expr.tumor), rownames(BRCA.Modules.OE)))

    # Define model input
    #TCGA.expr.tumor <- TCGA.expr.tumor[, match(TCGA.meta.tumor$sampleID, colnames(TCGA.expr.tumor))]

    Inputdata.centralGene <- t(TCGA.expr.tumor[rownames(TCGA.expr.tumor) %in% Common.CentralGenes, ]) %>% scale(center = TRUE, scale = TRUE) %>%
        as.data.frame()

    if(train_meta_data_type=="TCGA"){
        Inputdata.centralGene <- Inputdata.centralGene %>%
            mutate(Class = as.factor(TCGA.meta.tumor$SUBTYPE[match(rownames(Inputdata.centralGene), TCGA.meta.tumor$sampleID)]))
    }
    if(train_meta_data_type=="METABRIC"){
        Inputdata.centralGene <- Inputdata.centralGene %>%
            mutate(Class = as.factor(TCGA.meta.tumor$CLAUDIN_SUBTYPE[match(rownames(Inputdata.centralGene), TCGA.meta.tumor$PATIENT_ID)]))
    }


    Inputdata.centralGene <- Inputdata.centralGene[! is.na(Inputdata.centralGene$Class), ]

    CentralGenes.model <- fn_RF_classifier(Inputdata.centralGene, n.folds, repetitions, metric = Metric, tunegrid)

    #Test classification performance on second cohort
    Inputdata.centralGene.Metabric <- t(METABRIC.expr[Common.CentralGenes, ]) %>% scale(center = TRUE, scale = TRUE)
    Prediction.centralGenes.model <- predict(CentralGenes.model$Model, Inputdata.centralGene.Metabric)

    if(test_meta_data_type == "TCGA"){
        CentralGenes.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.centralGenes.model), as.factor(METABRIC.meta$SUBTYPE))
    }
    if(test_meta_data_type == "METABRIC"){
        CentralGenes.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.centralGenes.model), as.factor(METABRIC.meta$CLAUDIN_SUBTYPE))
    }

    return(CentralGenes.model)

}
#' build random classifiers
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param sponge_modules result of define_modules()
#' @param train_gene_expr expression data of train dataset,
#' genenames must be in rownames
#' @param test_gene_expr expression data of test dataset,
#' genenames must be in rownames
#' @param train_meta_data meta data of train dataset
#' @param test_meta_data meta data of test dataset
#' @param train_meta_data_type TCGA or METABRIC
#' @param test_meta_data_type TCGA or METABRIC
#' @param metric metric (Exact_match, Accuracy) (default: Exact_match)
#' @param tunegrid_c defines the grid for the hyperparameter optimization during
#' cross validation (caret package) (default: 1:100)
#' @param n.folds number of folds to be calculated
#' @param repetitions  number of k-fold cv iterations (default: 3)
#'
#' @param bin.size bin size (default: 100)
#' @param min.size minimum module size (default: 10)
#' @param max.size maximum module size (default: 200)
#' @param min.expr minimum expression (default: 10)
#' @param method Enrichment to be used (Overall Enrichment: OE or Gene Set
#' Variation Analysis: GSVA) (default: OE)
#' @param cores number of cores to be used to calculate entichment scores with gsva or ssgsea methods. Default 1
#' @param replace Possibility of keeping or removing (default) central
#' genes in the modules (default: F)
#'
#' @export
#'
#' @return randomized prediction model
# build_classifier_random<-function(sponge_modules,
#                                   train_gene_expr,
#                                   test_gene_expr,
#                                   train_meta_data,
#                                   test_meta_data,
#                                   train_meta_data_type="TCGA",
#                                   test_meta_data_type="TCGA",
#                                   metric="Exact_match",
#                                   tunegrid_c=c(1:100),
#                                   n.folds = 10,
#                                   repetitions=3,
#                                   min.size = 10,
#                                   bin.size = 100,
#                                   max.size = 200,
#                                   min.expression=10,
#                                   replace = F,
#                                   method = "OE"){
#
#     Sponge.modules<- sponge_modules
#     Metric<-metric
#
#     Sizes.modules <- lengths(Sponge.modules)
#
#     TCGA.expr.tumor<-train_gene_expr
#     METABRIC.expr<-test_gene_expr
#     TCGA.meta.tumor<-train_meta_data
#     METABRIC.meta<-test_meta_data
#
#     tunegrid<-expand.grid(.mtry=tunegrid_c)
#
#     # Define random modules
#     Random.Modules <- list()
#
#     for(j in 1:length(Size.modules)) {
#         Module.Elements <-  sample.int(n = nrow(TCGA.expr.tumor), size = Size.modules[j],
#                                        replace = replace)
#         # print(Module.Elements)
#         Random.Modules[[j]] <- rownames(TCGA.expr.tumor)[Module.Elements]
#     }
#     names(Random.Modules) <- names(Sponge.modules)
#
#     # Calculate enrichment scores for BRCA and METABRIC
#     BRCA.RandomModules.OE <- enrichment_modules(TCGA.expr.tumor, Random.Modules,
#                                                 bin.size = bin.size, min.size = min.size, max.size = max.size, method = method, min.expr =  min.expression)
#     METABRIC.RandomModules.OE <- enrichment_modules(METABRIC.expr, Random.Modules,
#                                                     bin.size = bin.size, min.size = min.size, max.size = max.size, method = method, min.expr =  min.expression)
#
#     #Find common modules
#
#     CommonModules <- intersect(rownames(BRCA.RandomModules.OE), rownames(METABRIC.RandomModules.OE))
#     BRCA.RandomModules.OE <- BRCA.RandomModules.OE[CommonModules, ]
#     METABRIC.RandomModules.OE <- METABRIC.RandomModules.OE[CommonModules, ]
#
#     # Train model
#     Inputdata.model.RANDOM <- t(BRCA.RandomModules.OE) %>% scale(center = T, scale = T) %>%
#         as.data.frame()
#
#     if(train_meta_data_type=="TCGA"){
#         Inputdata.model.RANDOM <- Inputdata.model.RANDOM %>%
#             mutate(Class = as.factor(TCGA.meta.tumor$SUBTYPE[match(rownames(Inputdata.model.RANDOM), TCGA.meta.tumor$sampleID)]))
#     }
#     if(train_meta_data_type=="METABRIC"){
#         Inputdata.model.RANDOM <- Inputdata.model.RANDOM %>%
#             mutate(Class = as.factor(TCGA.meta.tumor$CLAUDIN_SUBTYPE[match(rownames(Inputdata.model.RANDOM), TCGA.meta.tumor$PATIENT_ID)]))
#     }
#
#     Inputdata.model.RANDOM <- Inputdata.model.RANDOM[! is.na(Inputdata.model.RANDOM$Class), ]
#
#     Random.model <- fn_RF_classifier(Inputdata.model.RANDOM, n.folds, repetitions, metric = Metric, tunegrid)
#
#     #Test classification performance on second cohort
#     METABRIC.RandomModules.OE <- METABRIC.RandomModules.OE[ ,complete.cases(t(METABRIC.RandomModules.OE))]
#
#     if(test_meta_data_type=="TCGA"){
#         Meta.metabric <- METABRIC.meta$SUBTYPE[match(colnames(METABRIC.RandomModules.OE), METABRIC.meta$sampleID)]
#     }
#     if(test_meta_data_type=="METABRIC"){
#         Meta.metabric <- METABRIC.meta$CLAUDIN_SUBTYPE[match(colnames(METABRIC.RandomModules.OE), METABRIC.meta$PATIENT_ID)]
#     }
#
#     Input.Metabric.random <- t(METABRIC.RandomModules.OE) %>% scale(center = T, scale = T)
#     Prediction.random.model <- predict(Random.model$Model, Input.Metabric.random)
#     Random.model$Prediction.model <- Prediction.random.model
#     Random.model$ConfusionMatrix_testing <- confusionMatrix(as.factor(Prediction.random.model), Meta.metabric)
#
#     return(Random.model)
# }

#' Define random modules
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#'
#' @param sponge_modules result of define_modules()
#' @param train_gene_expr expression data of train dataset,
#' genenames must be in rownames
#' @param test_gene_expr expression data of test dataset,
#' genenames must be in rownames
#' @param train_meta_data meta data of train dataset
#' @param test_meta_data meta data of test dataset
#' @param train_meta_data_type TCGA or METABRIC
#' @param test_meta_data_type TCGA or METABRIC
#' @param metric metric (Exact_match, Accuracy) (default: Exact_match)
#' @param tunegrid_c defines the grid for the hyperparameter optimization during
#' cross validation (caret package) (default: 1:100)
#' @param n.folds number of folds to be calculated
#' @param repetitions  number of k-fold cv iterations (default: 3)
#'
#' @param bin.size bin size (default: 100)
#' @param min.size minimum module size (default: 10)
#' @param max.size maximum module size (default: 200)
#' @param min.expr minimum expression (default: 10)
#' @param method Enrichment to be used (Overall Enrichment: OE or Gene Set
#' Variation Analysis: GSVA) (default: OE)
#' @param cores number of cores to be used to calculate entichment scores with gsva or ssgsea methods. Default 1
#' @param replace Possibility of keeping or removing (default) central
#' genes in the modules (default: F)
#'
#' @export
#'
#' @return randomized prediction model
Random_spongEffects<-function(sponge_modules,
                              gene_expr,
                              min.size = 10,
                              bin.size = 100,
                              max.size = 200,
                              min.expression=10,
                              replace = F,
                              method = "OE",
                              cores = 1){

    Size.modules = sapply(sponge_modules, length)

    # Define random modules
    Random.Modules <- list()

    for(j in 1:length(Size.modules)) {
        Module.Elements <-  sample.int(n = nrow(gene_expr), size = Size.modules[j],
                                       replace = replace)
        # print(Module.Elements)
        Random.Modules[[j]] <- rownames(gene_expr)[Module.Elements]
    }
    names(Random.Modules) <- names(Sponge.modules)

    # Calculate enrichment scores with chosen method
    Random.modules <- enrichment_modules(gene_expr, Random.Modules,
                                         bin.size = bin.size, min.size = min.size, max.size = max.size, method = method, min.expr =  min.expression, cores)

    return(Random.modules)
}

#' plots the top x gini index modules (see Boniolo and Hoffmann 2022 et al. Figure 5)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import ComplexHeatmap
#' @import ggplot2
#' @import MetBrewer
#'
#' @param trained_model returned from train_and_test_model
#' @param k_modules top k modules to be shown (default: 25)
#' @param k_modules_red top k modules shown in red - NOTE: must be smaller
#' than k_modules (default: 10)
#' @param text_size text size (default 16)
#' @param bioMart_gene_symbol_columns bioMart dataset column for gene symbols
#' (e.g. human: hgnc_symbol, mouse: mgi_symbol)
#' (default: hgnc_symbol)
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#'
#' @export
#'
#' @return plot object for lollipop plot
plot_top_modules <- function(trained_model,
                             k_modules = 25,
                             k_modules_red = 10,
                             text_size = 16,
                             bioMart_gene_symbol_columns = "hgnc_symbol",
                             bioMart_gene_ensembl = "hsapiens_gene_ensembl") {

    SpongingActiivty.model <-trained.model$SpongingActivity.model
    trained.model<-trained_model

    final.model <- SpongingActiivty.model$Model$finalModel
    Variable.importance <- importance(final.model) %>% as.data.frame() %>%
        arrange(desc(MeanDecreaseGini)) %>%
        tibble::rownames_to_column('Module')

    mart <- useMart("ensembl", dataset=bioMart_gene_ensembl)
    genes <- Variable.importance$Module
    G_list <- getBM(filters = "ensembl_gene_id",
                    attributes = c("ensembl_gene_id",bioMart_gene_symbol_columns,'external_gene_name', "description"),
                    values = genes, mart = mart, useCache = FALSE)

    Variable.importance <- Variable.importance %>%
        mutate(Hugo_Symbol = G_list$hgnc_symbol[match(Variable.importance$Module, G_list$ensembl_gene_id)],
               Name_lncRNA_Results = ifelse(Hugo_Symbol == "", Module, Hugo_Symbol))

    grey_modules=k_modules-k_modules_red

    p<-Variable.importance[1:k_modules, ] %>%
        mutate(Analysed = c(rep("1", k_modules_red), rep("0",grey_modules))) %>%
        arrange(desc(MeanDecreaseGini)) %>%
        ggplot(aes(x=reorder(Name_lncRNA_Results, MeanDecreaseGini), y=MeanDecreaseGini)) +
        geom_point() +
        geom_segment( aes(x=Name_lncRNA_Results, xend=Name_lncRNA_Results, y=0, yend=MeanDecreaseGini, color = Analysed)) +
        scale_colour_manual(values=c("red", "black"), breaks = c("1", "0")) +
        coord_flip()+
        xlab("Module") +
        ylab("Mean decrease in Gini index") +
        theme_light()+
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_blank(),
            legend.position="none",
            legend.background = element_blank(),
            legend.direction="horizontal",
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            # axis.text=element_text(size=12),
            # axis.title=element_text(size=12), #,face="bold"
            # legend.text=element_text(size=12),
            text=element_text(size=16)
        )

    return(p)

}

#' plots the density of the model scores for subtypes (see Boniolo and Hoffmann 2022 et al. Fig. 2)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import ComplexHeatmap
#' @import ggplot2
#' @import MetBrewer
#' @import tidyr
#'
#' @param trained_model returned from train_and_test_model
#' @param modules output of enrichment_modules()
#' @param meta_data metadata of samples
#' (retrieved from prepare_tcga_for_spongEffects() or
#' from prepare_metabric_for_spongEffects())
#' @param meta_data_type TCGA or METABRIC
#' @param subtypes array of subtypes
#' (e.g., c("Normal", "LumA", "LumB", "Her2", "Basal"))
#' @param bioMart_gene_symbol_columns bioMart dataset column for gene symbols
#' (e.g. human: hgnc_symbol, mouse: mgi_symbol)
#' (default: hgnc_symbol)
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#'
#' @export
#'
#' @return plots density scores for subtypes
plot_density_scores <- function(trained_model,
                                modules,
                                meta_data,
                                meta_data_type="TCGA",
                                subtypes,
                                bioMart_gene_symbol_columns = "hgnc_symbol",
                                bioMart_gene_ensembl = "hsapiens_gene_ensembl"){
    ##FIGURE 2
    trained.model<-trained_model
    data_type<-meta_data_type

    title<-paste0("spongEffects - ",data_type," (Training)")

    SpongingActiivty.model <-trained.model$SpongingActivity.model

    final.model <- SpongingActiivty.model$Model$finalModel
    Variable.importance <- importance(final.model) %>% as.data.frame() %>%
        arrange(desc(MeanDecreaseGini)) %>%
        tibble::rownames_to_column('Module')

    mart <- useMart("ensembl", dataset=bioMart_gene_ensembl)
    genes <- Variable.importance$Module
    G_list <- getBM(filters = "ensembl_gene_id",
                    attributes = c("ensembl_gene_id",bioMart_gene_symbol_columns,'external_gene_name', "description"),
                    values = genes, mart = mart, useCache = FALSE)

    Variable.importance <- Variable.importance %>%
        mutate(Hugo_Symbol = G_list$hgnc_symbol[match(Variable.importance$Module, G_list$ensembl_gene_id)],
               Name_lncRNA_Results = ifelse(Hugo_Symbol == "", Module, Hugo_Symbol))

    if (data_type=="TCGA") {

        BRCA.Modules.OE<-modules
        TCGA.meta.tumor=meta_data

        TCGA.meta.tumor$SUBTYPE <- factor(TCGA.meta.tumor$SUBTYPE, levels = subtypes)
        DrivingModules.BRCA <- BRCA.Modules.OE[rownames(BRCA.Modules.OE) %in% Variable.importance$Module, ] %>%
            gather(Patient, Score) %>%
            mutate(Class = TCGA.meta.tumor$SUBTYPE[match(Patient, TCGA.meta.tumor$sampleID)]) %>%
            ggplot(aes(x = Score, y = Class, fill = Class)) +
            geom_density_ridges() +
            xlab(title) +
            scale_fill_manual(values=met.brewer("Hokusai2", 5)) +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = "none",
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  axis.text=element_text(size=25), axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = "black"),
                  axis.title.y =element_text(size=30),
                  axis.title.x =element_text(size=30))
        return(DrivingModules.BRCA)
    }

    if(data_type=="METABRIC"){

        METABRIC.Modules.OE<-modules
        METABRIC.meta=meta_data

        METABRIC.meta$CLAUDIN_SUBTYPE <- factor(METABRIC.meta$CLAUDIN_SUBTYPE, levels = subtypes)
        DrivingModules.METABRIC <- METABRIC.Modules.OE[rownames(METABRIC.Modules.OE) %in% Variable.importance$Module, ] %>%
            gather(Patient, Score) %>%
            mutate(Class = METABRIC.meta$CLAUDIN_SUBTYPE[match(Patient, METABRIC.meta$PATIENT_ID)]) %>%
            ggplot(aes(x = Score, y = Class, fill = Class)) +
            geom_density_ridges() +
            xlab(title) +
            scale_fill_manual(values=met.brewer("OKeeffe2", 5)) +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = "none",
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  axis.text=element_text(size=25), axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = "black"),
                  axis.title.y =element_text(size=30),
                  axis.title.x =element_text(size=30))

        return(DrivingModules.METABRIC)
    }




}
#' list of plots for (1) accuracy and (2) sensitivity + specificity
#' (see Boniolo and Hoffmann 2022 et al. Fig. 3a and Fig. 3b)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import ComplexHeatmap
#' @import ggplot2
#' @import MetBrewer
#'
#' @param trained_model returned from train_and_test_model
#' @param central_genes_model returned from build_classifier_central_genes()
#' @param random_model returned from train_and_test_model using the randomization
#' @param training_dataset_name name of training (e.g., TCGA)
#' @param testing_dataset_name name of testing set (e.g., METABRIC)
#' @param subtypes array of subtypes
#' (e.g., c("Normal", "LumA", "LumB", "Her2", "Basal"))
#'
#' @export
#'
#' @return list of plots for (1) accuracy and (2) sensitivity + specificity
plot_accuracy_sensitivity_specificity <- function(trained_model,
                                                  central_genes_model,
                                                  random_model,
                                                  training_dataset_name="TCGA",
                                                  testing_dataset_name="TCGA",
                                                  subtypes){
    trained.model<-trained_model
    CentralGenes.model<-central_genes_model
    Random.model<-random_model

    training_string<-paste0(training_dataset_name," (Training)")
    testing_string<-paste0(testing_dataset_name," (Testing)")
    set_names<-c(training_string,testing_string)

    SpongingActiivty.model <-trained.model$SpongingActivity.model
    #CentralGenes.model
    #Random.model

    ##FIGURE 3A
    Accuracy.df <- data.frame(Run = rep(set_names, 3),
                              Model = c(rep("Modules", 2), rep("Central Genes", 2), rep("Random", 2)),
                              Accuracy = c(SpongingActiivty.model$ConfusionMatrix_training[["overall"]][['Accuracy']], SpongingActiivty.model$ConfusionMatrix_testing[["overall"]][['Accuracy']],
                                           CentralGenes.model$ConfusionMatrix_training[["overall"]][['Accuracy']],CentralGenes.model$ConfusionMatrix_testing[["overall"]][['Accuracy']],
                                           Random.model$ConfusionMatrix_training[["overall"]][['Accuracy']], Random.model$ConfusionMatrix_testing[["overall"]][['Accuracy']]))

    Accuracy.df$Model <- factor(Accuracy.df$Model, levels = c("Modules", "Random", "Central Genes"))
    Accuracy.df$Run <- factor(Accuracy.df$Run, levels = set_names)

    Accuracy.plot <- Accuracy.df %>%
        ggplot(aes(x=Accuracy, y=Model)) +
        geom_line(aes(group = Model), color = c( "#17154f",  "#355828", "#bf3729",  "#17154f", "#355828", "#bf3729"), size = 1.5) + # Renoir, 3
        geom_point(aes(shape = Run), size = 6, color = c("#17154f", "#17154f",  "#355828","#355828", "#bf3729", "#bf3729")) +
        # scale_color_manual(values = c("#0a3351", "#b9563f", "#0a3351", "#b9563f","#0a3351", "#b9563f")) +
        scale_shape_manual(values=c(16, 17))+
        theme_light() +
        xlab("Subset Accuracy") +
        ylab("") +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_blank(),
            legend.position="top",
            legend.background = element_blank(),
            legend.direction="horizontal",
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.text=element_text(size=15),
            axis.title=element_text(size=20),
            legend.text=element_text(size=12), axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = c( "#17154f",  "#bf3729", "#355828"))
        )

    # Sensitivity and Specificity
    # Figure 3b
    Metrics.SpongeModules.training <- SpongingActiivty.model$ConfusionMatrix_training[["byClass"]][c(1:5), c(1,2)] %>% as.data.frame() %>%
        mutate(Model = "Modules") %>% tibble::rownames_to_column('Class') %>% gather(Metric, Value, Sensitivity:Specificity,-Class, -Model)
    Metrics.Random.training <- Random.model$ConfusionMatrix_training[["byClass"]][c(1:5), c(1,2)] %>% as.data.frame() %>%
        mutate(Model = "Random") %>% tibble::rownames_to_column('Class') %>% gather(Metric, Value, Sensitivity:Specificity,-Class, -Model)
    Metrics.CentralGenes.training <- CentralGenes.model$ConfusionMatrix_training[["byClass"]][c(1:5), c(1,2)] %>% as.data.frame()%>%
        mutate(Model = "Central Genes") %>% tibble::rownames_to_column('Class') %>% gather(Metric, Value, Sensitivity:Specificity,-Class, -Model)

    Metrics.training <- rbind(Metrics.SpongeModules.training, rbind(Metrics.Random.training, Metrics.CentralGenes.training)) %>%
        mutate(Run = training_string)

    Metrics.SpongeModules.testing <- SpongingActiivty.model$ConfusionMatrix_testing[["byClass"]][c(1:5), c(1,2)] %>% as.data.frame() %>%
        mutate(Model = "Modules") %>% tibble::rownames_to_column('Class') %>% gather(Metric, Value, Sensitivity:Specificity,-Class, -Model)
    Metrics.Random.testing <- Random.model$ConfusionMatrix_testing[["byClass"]][c(1:5), c(1,2)] %>% as.data.frame()%>%
        mutate(Model = "Random") %>% tibble::rownames_to_column('Class') %>% gather(Metric, Value, Sensitivity:Specificity,-Class, -Model)
    Metrics.CentralGenes.testing <- CentralGenes.model$ConfusionMatrix_testing[["byClass"]][c(1:5), c(1,2)] %>% as.data.frame()%>%
        mutate(Model = "Central Genes") %>% tibble::rownames_to_column('Class') %>% gather(Metric, Value, Sensitivity:Specificity,-Class, -Model)

    Metrics.testing <- rbind(Metrics.SpongeModules.testing, rbind(Metrics.Random.testing, Metrics.CentralGenes.testing)) %>%
        mutate(Run = testing_string)

    Metrics <- rbind(Metrics.training, Metrics.testing)
    Metrics$Model <- factor(Metrics$Model, levels = c("Modules", "Random", "Central Genes"))
    Metrics$Run <- factor( Metrics$Run, levels = set_names)

    Metrics$Class <- gsub("Class: ", "", Metrics$Class)
    Metrics$Class <- factor(Metrics$Class, levels =subtypes)

    Metrics.plot <- Metrics %>%
        ggplot(aes(x = Metric, y = Value, fill = Model)) +
        geom_bar(position = "dodge", stat = "identity", width = 0.5) +
        facet_grid(Run ~ Class) +
        scale_fill_manual(values=c( "#17154f",  "#bf3729", "#355828"),
                          breaks = c("Modules", "Random", "Central Genes")) +
        xlab("Metric") +
        ylab("") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill="white", size=1, color="black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              strip.text.x = element_text(size=15,  color="black"),
              strip.text.y = element_text(size=15, color="black"),
              axis.text=element_text(size=15),
              axis.title=element_text(size=20),
              legend.text=element_text(size=15), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, color = "black"),
              axis.text.y = element_text(color = "black"),
              legend.title=element_text(size=15))

    metric_plots<-list(Accuracy.plot,Metrics.plot)
    names(metric_plots) <- c("Accuracy.plot","Metrics.plot")

    return(metric_plots)
}

#' plots the confusion matrix from spongEffects train_and_test()
#' (see Boniolo and Hoffmann 2022 et al. Fig. 3a and Fig. 3b)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import ComplexHeatmap
#' @import ggplot2
#' @import MetBrewer
#'
#' @param trained_model returned from train_and_test_model
#' @param subtypes_testing_factors subtypes of testing samples as factors
#' @return plot of the confusion matrix
#'
#' @export
#'
#' @return returns confusion matrix plots of the trained model
plot_confusion_matrices <- function(trained_model,
                                    subtypes.testing.factors){

    trained.model<-trained_model
    # Visualize confusion matrices  -------------------------------------------
    # SpongEffect model (Supplementary Confusion Matrices)

    Meta.metabric<- na.omit(subtypes.testing.factors)

    # Prediction.SpongeModules<- as.character(trained.model$Prediction.model)

    Prediction.SpongeModules <- data.frame(Predicted = as.factor(trained.model$Prediction.model),
                                           Observed = as.factor(subtypes.testing.factors))
    Prediction.SpongeModules.cm <- confusion_matrix(targets = Prediction.SpongeModules$Observed,
                                                    predictions = Prediction.SpongeModules$Predicted)

    Confusion.matrix.SpongeModules <- cvms::plot_confusion_matrix(Prediction.SpongeModules.cm$`Confusion Matrix`[[1]],
                                                                  add_sums = T,
                                                                  add_normalized = FALSE,
                                                                  add_col_percentages = FALSE,
                                                                  add_row_percentages = FALSE,
                                                                  sums_settings = sum_tile_settings(
                                                                      palette = "Oranges",
                                                                      label = "Total",
                                                                      tc_tile_border_color = "black"),
                                                                  darkness = 0)

    return(Confusion.matrix.SpongeModules)
}

#' plots the heatmaps from training_and_test_model
#' (see Boniolo and Hoffmann 2022 et al. Fig. 6)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import ComplexHeatmap
#' @import ggplot2
#' @import MetBrewer
#'
#' @param trained_model returned from train_and_test_model
#' @param sponge_modules result of define_modules()
#' @param enrichment_modules return of enrichment_modules()
#' @param meta_data meta data
#' @param meta_data_type TCGA or METABRIC
#' @param subtypes array of subtypes
#' (e.g., c("Normal", "LumA", "LumB", "Her2", "Basal"))
#' @param bioMart_gene_symbol_columns bioMart dataset column for gene symbols
#' (e.g. human: hgnc_symbol, mouse: mgi_symbol)
#' (default: hgnc_symbol)
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#' @param met_brewer.palette_subtypes which met brewer colour palette should be
#' used for the subtypes (default: "Austria")
#' @param met_brewer.palette_heatmap which  met brewer colour palette should be
#' used for the heatmap (default "Troy")
#'
#' @export
#'
#' @return
#' NOT FUNCTIONAL
plot_heatmaps_training_test<-function(trained_model,
                                      sponge_modules,
                                      enrichment_modules,
                                      meta_data,
                                      meta_data_type = "TCGA",
                                      subtypes,
                                      bioMart_gene_symbol_columns = "hgnc_symbol",
                                      bioMart_gene_ensembl = "hsapiens_gene_ensembl",
                                      met_brewer.palette_subtypes = "Austria",
                                      met_brewer.palette_heatmap = "Troy"){

    trained.model<-trained_model
    Sponge.modules<-sponge_modules
    enrichment.modules<-enrichment_modules
    meta.data<-meta_data
    type.data.set<-meta_data_type

    colors_palette <- met.brewer("Austria",n=length(subtypes))
    colors_palette <- as.character(colors_palette)
    names(colors_palette)<-subtypes
    Colours.subtypes <- colors_palette

    METABRIC.meta<-meta.data
    METABRIC.Modules.OE<-enrichment.modules

    # Variable importance -----------------------------------------------------
    SpongingActiivty.model <-trained.model$SpongingActivity.model

    final.model <- SpongingActiivty.model$Model$finalModel
    Variable.importance <- importance(final.model) %>% as.data.frame() %>%
        arrange(desc(MeanDecreaseGini)) %>%
        tibble::rownames_to_column('Module')

    mart <- useMart("ensembl", dataset=bioMart_gene_ensembl)
    genes <- Variable.importance$Module
    G_list <- getBM(filters = "ensembl_gene_id",
                    attributes = c("ensembl_gene_id",bioMart_gene_symbol_columns,'external_gene_name', "description"),
                    values = genes, mart = mart, useCache = FALSE)

    Variable.importance <- Variable.importance %>%
        mutate(Hugo_Symbol = G_list$hgnc_symbol[match(Variable.importance$Module, G_list$ensembl_gene_id)],
               Name_lncRNA_Results = ifelse(Hugo_Symbol == "", Module, Hugo_Symbol))

    # rownames(BRCA.Modules.OE) <- Variable.importance$Name_lncRNA_Results[match(rownames(BRCA.Modules.OE), Variable.importance$Module)]
    # Visualize Heatmaps  ---------------------------------------------------------
    # Define colours heatmaps

    ##FIGURE 6

    col <- rev(met.brewer(met_brewer.palette_heatmap,n=100,type="continuous"))

    # METABRIC
    # Define annotation layer

    if(type.data.set=="METABRIC"){
        METABRIC.meta.p <- METABRIC.meta[match(colnames(METABRIC.Modules.OE), METABRIC.meta$PATIENT_ID), ]
        METABRIC.meta.p$CLAUDIN_SUBTYPE <- factor(METABRIC.meta.p$CLAUDIN_SUBTYPE, levels = subtypes)
        HDHAnnotation <- HeatmapAnnotation(Subtype = METABRIC.meta.p$CLAUDIN_SUBTYPE,
                                           annotation_name_gp = 15,
                                           # Er_IHC = METABRIC.meta.p$ER_IHC,
                                           # Her_SNP = METABRIC.meta.p$HER2_SNP6,
                                           col = list(Subtype = Colours.subtypes))

        Heatmap.METABRIC <- METABRIC.Modules.OE[match(Variable.importance$Module[1:length(Variable.importance$Module)], rownames(METABRIC.Modules.OE)), ]
        rownames(Heatmap.METABRIC) <- Variable.importance$Name_lncRNA_Results[1:length(Variable.importance$Module)]

        # Plot
        Heatmap.METABRIC.p <- Heatmap.METABRIC %>%
            t() %>% scale() %>% t() %>%
            Heatmap(show_row_names = T, show_column_names = F,  top_annotation = HDHAnnotation, column_split = METABRIC.meta.p$CLAUDIN_SUBTYPE,
                    cluster_column_slices = F,
                    column_title_gp = gpar(fontsize = 20),
                    row_names_gp = gpar(fontsize = 15),
                    row_title_gp = gpar(fontsize = 20),
                    row_km = 3,
                    col = col,
                    heatmap_legend_param = list(
                        title = "Z score", at = c(-3, 0, 3),
                        labels = c("-3", "0", "3"), direction = "horizontal", legend_width = unit(3, "cm"),
                        labels_gp = gpar(fontsize = 15),
                        title_gp = gpar(fontsize = 15),
                        title_position = "topcenter"
                    ),  show_parent_dend_line = FALSE,
                    # row_title = "Modules",
                    column_title = type.data.set,width = unit(16, "cm"), height = unit(16, "cm"), show_heatmap_legend = FALSE)
        return(Heatmap.METABRIC.p)
    }

    if(type.data.set=="TCGA"){

        METABRIC.meta.p <- METABRIC.meta[match(colnames(METABRIC.Modules.OE), METABRIC.meta$sampleID), ]
        METABRIC.meta.p$SUBTYPE <- factor(METABRIC.meta.p$SUBTYPE, levels = subtypes)
        HDHAnnotation <- HeatmapAnnotation(Subtype = METABRIC.meta.p$SUBTYPE,
                                           annotation_name_gp = 15,
                                           # Er_IHC = METABRIC.meta.p$ER_IHC,
                                           # Her_SNP = METABRIC.meta.p$HER2_SNP6,
                                           col = list(Subtype = Colours.subtypes))

        Heatmap.METABRIC <- METABRIC.Modules.OE[match(Variable.importance$Module[1:length(Variable.importance$Module)], rownames(METABRIC.Modules.OE)), ]
        rownames(Heatmap.METABRIC) <- Variable.importance$Name_lncRNA_Results[1:length(Variable.importance$Module)]

        # Plot
        Heatmap.METABRIC.p <- Heatmap.METABRIC %>%
            t() %>% scale() %>% t() %>%
            Heatmap(show_row_names = T, show_column_names = F,  top_annotation = HDHAnnotation, column_split = METABRIC.meta.p$SUBTYPE,
                    cluster_column_slices = F,
                    column_title_gp = gpar(fontsize = 20),
                    row_names_gp = gpar(fontsize = 15),
                    row_title_gp = gpar(fontsize = 20),
                    row_km = 3,
                    col = col,
                    heatmap_legend_param = list(
                        title = "Z score", at = c(-3, 0, 3),
                        labels = c("-3", "0", "3"), direction = "horizontal", legend_width = unit(3, "cm"),
                        labels_gp = gpar(fontsize = 15),
                        title_gp = gpar(fontsize = 15),
                        title_position = "topcenter"
                    ),  show_parent_dend_line = FALSE,
                    # row_title = "Modules",
                    column_title = type.data.set,width = unit(16, "cm"), height = unit(16, "cm"), show_heatmap_legend = FALSE)
        return(Heatmap.METABRIC.p)
    }
}

#' plots the heatmap of miRNAs invovled in the interactions of the modules
#' (see Boniolo and Hoffmann 2022 et al. Fig. 7a)
#'
#' @import tidyverse
#' @import caret
#' @import dplyr
#' @import Biobase
#' @import biomaRt
#' @import randomForest
#' @import ggridges
#' @import cvms
#' @import miRBaseConverter
#' @import ComplexHeatmap
#' @import ggplot2
#' @import MetBrewer
#' @import grid
#'
#' @param sponge_modules result of define_modules()
#' @param trained_model returned from train_and_test_model
#' @param gene_mirna_candidates output of SPONGE or SPONGEdb (miRNAs_significance)
#' @param k_modules top k modules to be shown (default: 25)
#' @param filter_miRNAs min rowsum to be reach of miRNAs (default: 3.0)
#' @param bioMart_gene_symbol_columns bioMart dataset column for gene symbols
#' (e.g. human: hgnc_symbol, mouse: mgi_symbol)
#' (default: hgnc_symbol)
#' @param bioMart_gene_ensembl bioMart gene ensemble name
#' (e.g., hsapiens_gene_ensembl).
#'
#' @export
#'
#' @return plot object
plot_involved_miRNAs_to_modules<-function(sponge_modules,
                                          trained_model,
                                          gene_mirna_candidates,
                                          k_modules = 25,
                                          filter_miRNAs = 3.0,
                                          bioMart_gene_symbol_columns = "hgnc_symbol",
                                          bioMart_gene_ensembl = "hsapiens_gene_ensembl"){

    Sponge.modules<-sponge_modules
    trained.model<-trained_model
    miRNAs_significance<-gene_mirna_candidates

    miRNAs_significance.downstream<-miRNAs_significance
    miRNAs_significance.downstream.Target<-miRNAs_significance

    SpongingActiivty.model <-trained.model$SpongingActivity.model

    final.model <- SpongingActiivty.model$Model$finalModel
    Variable.importance <- importance(final.model) %>% as.data.frame() %>%
        arrange(desc(MeanDecreaseGini)) %>%
        tibble::rownames_to_column('Module')

    Variable.importance.top_k=Variable.importance[1:k_modules, ]$Module
    Sponge.interesting.modules<-Sponge.modules[Variable.importance.top_k]
    Sponge.modules.downstrean<-Sponge.interesting.modules

    httr::set_config(httr::config(ssl_verifypeer = FALSE))
    mart <- useDataset(bioMart_gene_ensembl, useMart("ensembl"))

    vec_colnames = vector(length = k_modules)
    vec_colnames_ensg = names(Sponge.modules.downstrean)
    df_centralnodes_map = data.frame()
    df_mirnas_map = data.frame()

    count=0
    for (ensg in names(Sponge.modules.downstrean)) {
        count=count+1

        df_intern = data.frame(matrix(ncol=1,nrow=1))
        colnames(df_intern)<-c("Geneid")
        df_intern <- rbind(df_intern,ensg)
        df_intern <- df_intern[-c(1),]
        df_intern <- data.frame(df_intern)
        colnames(df_intern)<-c("Geneid")
        df_intern <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",bioMart_gene_symbol_columns,"description"),values=df_intern$Geneid,mart= mart)

        name <- df_intern$hgnc_symbol
        if(is.na(name))
        {
            name <- df_intern$ensembl_gene_id
        }
        vec_colnames[count] = name

        df_centralnodes_map<-rbind(df_centralnodes_map,df_intern)
    }

    df_heatmap = data.frame(matrix(ncol=k_modules,nrow=0))
    df_heatmap = as.data.frame(df_heatmap)
    colnames(df_heatmap)<-vec_colnames
    #length(names(df_heatmap))

    for (ensg in vec_colnames_ensg)
    {
        #ensg = "ENSG00000237233"
        print(ensg)
        module_symbole = ensg

        df_symbolname = df_centralnodes_map[which(df_centralnodes_map$ensembl_gene_id == ensg), ]

        module_symbole = df_symbolname$hgnc_symbol

        if(is.na(module_symbole))
        {
            module_symbole=ensg
        }

        vec_all_targets = Sponge.modules.downstrean[[ensg]]

        df_intern = data.frame(matrix(ncol=1,nrow=length(vec_all_targets)))
        colnames(df_intern)<-c("Geneid")
        df_intern$Geneid = vec_all_targets
        df_intern <- df_intern[-c(1),]
        df_intern <- data.frame(df_intern)
        colnames(df_intern)<-c("Geneid")
        df_intern <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",bioMart_gene_symbol_columns,"description"),values=df_intern$Geneid,mart= mart)


        df_mirnas_intern=data.frame(miRNAs_significance.downstream[[ensg]])
        df_mirnas_intern=miRNA_AccessionToName(df_mirnas_intern$mirna)

        df_mirna_counts = data.frame(matrix(ncol = length(df_mirnas_intern$Accession)))
        colnames(df_mirna_counts) <- df_mirnas_intern$TargetName
        df_mirna_counts[is.na(df_mirna_counts)] <- 0

        count_targets = 0

        df_mirnas_map<-rbind(df_mirnas_map,df_mirnas_intern)

        for (ensg_target in df_intern$ensembl_gene_id)
        {
            #ensg_target="ENSG00000237970"

            df_symbolname = df_intern[which(df_intern$ensembl_gene_id == ensg_target), ]

            #symbolname = df_symbolname$hgnc_symbol
            symbolname=ensg_target

            if(!is.na(symbolname))
            {

                miRNAs_thistarget <- miRNAs_significance.downstream.Target[[symbolname]]$mirna
                miRNAs_thistarget=miRNA_AccessionToName(miRNAs_thistarget)

                intersect_mirnas = intersect(miRNAs_thistarget$TargetName,df_mirnas_intern$TargetName)

                for (i_mirnas in intersect_mirnas)
                {
                    if(!is.na(i_mirnas)){
                        already_counted = df_mirna_counts[,i_mirnas]
                        already_counted = already_counted + 1
                        df_mirna_counts[,i_mirnas] = already_counted
                    }

                }
            }

            count_targets=count_targets+1
        }

        for (mirna_name in colnames(df_mirna_counts)) {
            #mirna_name="hsa-miR-4999-5p"
            if(!is.na(mirna_name)){
                if(!any(row.names(df_heatmap) == mirna_name))
                {
                    df_heatmap[nrow(df_heatmap)+1,] <- 0
                    rownames(df_heatmap)[length(rownames(df_heatmap))]=mirna_name
                }

                freq = df_mirna_counts[,mirna_name]/count_targets

                df_heatmap[mirna_name, module_symbole] = freq
            }

        }
    }

    df_heatmap<-df_heatmap[rowSums(df_heatmap[])>filter_miRNAs,]
    df_heatmap<-as.matrix(df_heatmap)

    col.heatmap <- met.brewer("VanGogh3",n=5,type="continuous")

    # Heatmap
    heatmap.miRNA <- df_heatmap %>%
        Heatmap(na_col = "white", col =  col.heatmap, # right_annotation = row.Annotation,
                show_row_names = T, show_column_names = T,
                heatmap_legend_param = list(
                    title = "Frequency",
                    at = c(0, 0.5, 1),
                    legend_height = unit(1.5, "cm"),
                    labels_gp = gpar(fontsize = 14),
                    title_gp = gpar(fontsize = 14),
                    direction = "horizontal", legend_width = unit(3, "cm"),
                    title_position = "topleft"),
                #top_annotation = column.Annotation,
                column_title = "Module", row_title = "miRNA",
                row_title_gp = gpar(fontsize = 20),
                column_title_gp = gpar(fontsize = 20),
                row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 15), column_names_rot = 45,
                rect_gp = gpar(col = "white", lwd = 0.5),
                width = unit(20, "cm"), height = unit(28, "cm"))

    return(heatmap.miRNA)

    #draw(heatmap.miRNA, heatmap_legend_side = "bottom",
    #     annotation_legend_side = "bottom")
}

