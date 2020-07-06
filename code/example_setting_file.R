
# to have compatible versions of Rdata
options(save.defaults = list(version = 2))  # for compatibility reasons


# in this file, settings that are specific for a run on a dataset

# gives path to output folder
pipOutFold <- <PATH_TO_OUTPUT_FOLDER>

# full path (starting with /mnt/...)
# following format expected for the input
# colnames = samplesID
# rownames = geneID
# !!! geneID are expected not difficulted

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 0_prepGeneData
# *************************************************************************************************************************

# UPDATE 07.12.2018: for RSEM data, the "analog" FPKM file is provided separately (built in prepData)
rna_fpkmDT_file <- <PATH TO FPKM EXPRESSION DATA>

rnaseqDT_file <- <PATH TO GENE EXPRESSION DATA>
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

removeDupGeneID <- TRUE # not tested with FALSE

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# labels for conditions
cond1 <- <NAME_CONDITION_1>
cond2 <- <NAME_CONDITION_2>

# path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
sample1_file <- <PATH_TO_RDATA_FILE_WITH_COND1_SAMPLES>
sample2_file <- <PATH_TO_RDATA_FILE_WITH_COND2_SAMPLES>

inputDataType <- "RSEM"

nCpu <- <NBR_AVAILABLE_CPU> # we used 25


# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
TADpos_file <- <PATH_TO_TAD_POSITIONS_FILE>
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

gene2tadDT_file <- <PATH_TO_TAD_GENE_TO_TAD_ASSIGNMENT_FILE>
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1


# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
# *************************************************************************************************************************


nRandomPermut <- <NBR_PERMUT> # we used 100000

min_counts <- <MIN_COUNT># we used 5
min_sampleRatio <- <MIN_RATIO_SAMPLES_THAT_SHOULD_HAVE_min_counts> # we used 0.8


