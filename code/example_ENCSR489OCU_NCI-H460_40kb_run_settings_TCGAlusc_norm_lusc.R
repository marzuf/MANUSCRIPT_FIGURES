# in this file, settings that are specific for a run on a dataset

# gives path to output folder
pipOutFold <- "EXAMPLE/OUTPUT_FOLDER/TCGAluad_norm_luad"

# full path (starting with /mnt/...)
# following format expected for the input
# colnames = samplesID
# rownames = geneID
# !!! geneID are expected not difficulted

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 0_prepGeneData
# *************************************************************************************************************************

# UPDATE 07.12.2018: for RSEM data, the "analog" FPKM file is provided separately (built in prepData)
# cp /mnt/ed4/marie/other_datasets/TCGAluad_norm_luad/fpkmDT.Rdata EXAMPLE/DATA/TCGAluad_norm_luad_fpkmDT.Rdata
rna_fpkmDT_file <- "EXAMPLE/DATA/TCGAluad_norm_luad_fpkmDT.Rdata"
# cp /mnt/ed4/marie/other_datasets/TCGAluad_norm_luad/rnaseqDT_v2.Rdata EXAMPLE/DATA/TCGAluad_norm_luad_rnaseqDT_v2.Rdata 
rnaseqDT_file <- "EXAMPLE/DATA/TCGAluad_norm_luad_rnaseqDT_v2.Rdata"
my_sep <- "\t"
# input is Rdata or txt file ?
# TRUE if the input is Rdata
inRdata <- TRUE

# can be ensemblID, entrezID, geneSymbol
geneID_format <- "entrezID"
stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

# are geneID rownames ? -> "rn" or numeric giving the column
geneID_loc <- "rn"
stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

removeDupGeneID <- TRUE

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# labels for conditions
cond1 <- "norm"
cond2 <- "luad"

# path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
# cp /mnt/ed4/marie/other_datasets/TCGAluad_norm_luad/norm_ID.Rdata EXAMPLE/DATA/norm_ID.Rdata
# cp /mnt/ed4/marie/other_datasets/TCGAluad_norm_luad/luad_ID.Rdata EXAMPLE/DATA/luad_ID.Rdata
sample1_file <- "EXAMPLE/DATA/norm_ID.Rdata"
sample2_file <- "EXAMPLE/DATA/luad_ID.Rdata"


inputDataType <- "RSEM"


# path to output folder:

# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
# cp /mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/ENCSR489OCU_NCI-H460_40kb/genes2tad/all_assigned_regions.txt EXAMPLE/DATA/ENCSR489OCU_NCI-H460_all_assigned_regions.txt
TADpos_file <- file.path("EXAMPLE/DATA/ENCSR489OCU_NCI-H460_all_assigned_regions.txt")
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

# cp /mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/ENCSR489OCU_NCI-H460_40kb/genes2tad/all_genes_positions.txt EXAMPLE/DATA/ENCSR489OCU_NCI-H460_all_genes_positions.txt
gene2tadDT_file <- file.path("EXAMPLE/DATA/ENCSR489OCU_NCI-H460_all_genes_positions.txt")
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1

# overwrite main_settings.R: nCpu <- 25
nCpu <- 40

# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
# *************************************************************************************************************************

# number of permutations
nRandomPermut <- 100000


min_counts <- 5
min_sampleRatio <- 0.8


# to have compatible versions of Rdata
options(save.defaults = list(version = 2))


# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATION DATA FOR CORRELATION
# *************************************************************************************************************************

all_permutCorr_data <- "/mnt/etemp/marie/PIPELINE/OUTPUT_FOLDER"
nbrCorrPermutCheck <- 58 
corrDiscardPattern <- "RANDOM|PERMUT"
refineMatchPattern <- "7sameNbr_"

#all_permutCorr_data <- "data/all_sample_corrValues.RData"
#nbrCorrPermutCheck <- 58 









