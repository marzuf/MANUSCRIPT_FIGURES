cat("load main_settings.R\n")

########################################################
### set paths to main files needed accross the workflow + general settings
########################################################

historyDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/data/gene_history_reordered.txt")

# file with mapping from entrez to chromosomic positions
#entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2position/filter_entrez_map.Rdata") # previous also held symbols
# holds only netrezID and positions
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")

# file with symbol synonyms to entrezID
#synoDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/all_entrez_syno_1.Rdata")
# mapping entrezID and all possible symbols
symbolDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_SYMBOL/final_entrez2syno.txt")

# mapping entrezID and ensemblID
ensemblDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_ENSEMBL/final_entrez2ensembl.txt")


nCpu <- 30

########################################################
### which kind of data to prepare and use across pipeline
########################################################

# if TRUE, from scripts 2_ to further, use only the genes that pass the DE analysis threshold for all analyses
useFilterCountData <- TRUE

# use only TADs (and corresponding genes)>= minSize et < maxQuantile and genes belonging to those TADs
useFilterSizeData <- TRUE

# if TRUE, from scripts 2_ to further, only focus on  TAD regions
useTADonly <- TRUE


########################################################
### TAD filter number of genes - for filter data (above) and used in 2_runWilcoxonTAD etc.
########################################################

# keep TAD with >= minNbrGeneTAD
minNbrGeneTAD <- 3

# keep TAD with <= maxQuant Genes !!! UPDATE HERE: FOR RUNNING PIPELINE V2 -> FOR TOPDOM, USE 0.99
maxQuantGeneTAD <- 0.99


########################################################
### permutation settings - 5_runPermutationsMedian.R        => STEP 5
########################################################

withExprClass <- TRUE
# number of classes of expression used for permutate expression data
permutExprClass <- 5


########################################################
### which quantile of the permutations to consider          => STEPS 8c, 12c, 14c, 17c
########################################################

permThresh <- 0.95


