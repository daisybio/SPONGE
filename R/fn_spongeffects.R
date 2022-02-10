

prepare_tcga_for_spongEffects
prepare_metabric_for_spongEffects
filter_ceRNA_network
get_lncRNA_modules
Define_Modules
# Module size distribution
Size.modules <- sapply(Sponge.modules, length)
BRCA.Modules.OE <- Enrichment_Modules(TCGA.expr.tumor, Sponge.modules,
                                      bin.size = 100, min.size = 10, max.size = 200, method = "OE")