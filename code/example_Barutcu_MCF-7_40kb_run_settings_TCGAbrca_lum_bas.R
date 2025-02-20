    # > file written: Sun, 09 Dec 2018 07:36:20 +0100 

    # in this file, settings that are specific for a run on a dataset

    # gives path to output folder
    pipOutFold <- "OUTPUT_FOLDER/TCGAbrca_lum_bas"

    # full path (starting with /mnt/...)
    # following format expected for the input
    # colnames = samplesID
    # rownames = geneID
    # !!! geneID are expected not difficulted

    # *************************************************************************************************************************
    # ************************************ SETTINGS FOR 0_prepGeneData
    # *************************************************************************************************************************

    # UPDATE 07.12.2018: for RSEM data, the "analog" FPKM file is provided separately (built in prepData)
    rna_fpkmDT_file <- "/mnt/ed4/marie/other_datasets/TCGAbrca_lum_bas/fpkmDT.Rdata"

    rnaseqDT_file <- "/mnt/ed4/marie/other_datasets/TCGAbrca_lum_bas/rnaseqDT_v2.Rdata"
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
    cond1 <- "lum"
    cond2 <- "bas"

    # path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
    sample1_file <- "/mnt/ed4/marie/other_datasets/TCGAbrca_lum_bas/lum_ID.Rdata"
    sample2_file <- "/mnt/ed4/marie/other_datasets/TCGAbrca_lum_bas/bas_ID.Rdata"

    inputDataType <- "RSEM"

    nCpu <- 20


# > file edited: Wed, 25 Mar 2020 18:10:30 +0100 

# path to output folder:
pipOutFold <- "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/Barutcu_MCF-7_40kb/TCGAbrca_lum_bas"

# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
TADpos_file <- paste0(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/Barutcu_MCF-7_40kb/genes2tad/all_assigned_regions.txt")
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

gene2tadDT_file <- paste0(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/Barutcu_MCF-7_40kb/genes2tad/all_genes_positions.txt")
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1

# overwrite main_settings.R: nCpu <- 25
nCpu <- 40

# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
# *************************************************************************************************************************


nRandomPermut <- 100000

# added here 13.08.2019 to change the filtering of min. read counts
if(exists("minCpmRatio")) rm(minCpmRatio)
min_counts <- 5
min_sampleRatio <- 0.8


# to have compatible versions of Rdata
options(save.defaults = list(version = 2))

