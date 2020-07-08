---
title: "Input data"
output: 
  html_document:
    self_contained: true
    keep_md: true
    toc: true
    toc_float: true
---



### Data with gene ID information

<ul>
<li>text files expected</li>
<li> see below for the expected format</li>
</ul>

##### Gene position information 


```r
# from main_settings.R
entrezDT_file <- paste0("gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
head(gff_dt)
#>    entrezID chromo start   end   assembly strand    symbol
#> 1 100287102   chr1 11874 14409 GRCh37.p13      +   DDX11L1
#> 2    653635   chr1 14362 29370 GRCh37.p13      -    WASH7P
#> 3 100302278   chr1 30366 30503 GRCh37.p13      + MIR1302-2
#> 4    645520   chr1 34611 36081 GRCh37.p13      -   FAM138A
#> 5     79504   chr1 52453 53396 GRCh37.p13      +    OR4G4P
#> 6    403263   chr1 63016 63885 GRCh37.p13      +   OR4G11P
```

##### Older entrez IDs mapping


```r
# from main_settings.R
historyDT_file <- file.path("gene_history_reordered.txt")
historyDT <- read.delim(historyDT_file,
                          stringsAsFactors = FALSE,
                          header=TRUE)
head(historyDT)
#>   entrezID mappingID
#> 1     7003         8
#> 2    51596        42
#> 3     1261        44
#> 4    54714        45
#> 5     7291        57
#> 6   643309        61
```

##### Gene symbol - entrez ID mapping


```r
# from main_settings.R
symbolDT_file <- paste0("final_entrez2syno.txt")
symbolDT <- read.delim(symbolDT_file,
                          stringsAsFactors = FALSE,
                          header=TRUE)
head(symbolDT)
#>   entrezID   symbol
#> 1        1     A1BG
#> 2        1      A1B
#> 3        1      ABG
#> 4        1      GAB
#> 5        1 HYST2477
#> 6        2      A2M
```

##### Gene ensemble ID - entrez ID mapping


```r
# from main_settings.R
ensemblDT_file <- paste0("final_entrez2ensembl.txt")
ensemblDT <- read.delim(ensemblDT_file,
                          stringsAsFactors = FALSE,
                          header=TRUE)
head(ensemblDT)
#>   entrezID       ensemblID
#> 1        1 ENSG00000121410
#> 2        2 ENSG00000175899
#> 3        3 ENSG00000256069
#> 4        9 ENSG00000171428
#> 5       10 ENSG00000156006
#> 6       12 ENSG00000196136
```

### TAD positions and gene-to-TAD assignment

<ul>
<li>text files expected</li>
<li>4-column TAD positions: chromosome/region/start/end</li>
<li>5-column gene information: geneID/chromosome/start/end/region</li>
<li>region label is expected to be like <em>chr<Nbr>_TAD<nbr></em></li>
</ul>

(order of the columns and extra-columns do not matter)


```r
# TAD positions
TADpos_file <- file.path("..", "EXAMPLE", "DATA", "ENCSR489OCU_NCI-H460_all_assigned_regions.txt")
tad_pos_dt <- read.delim(TADpos_file,
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("chromosome", "region", "start", "end"))
head(tad_pos_dt)
#>   chromosome       region   start     end
#> 1      chr10 chr10_BOUND1       1   80000
#> 2      chr10   chr10_TAD1   80001  720000
#> 3      chr10   chr10_TAD2  720001 1120000
#> 4      chr10   chr10_TAD3 1120001 1320000
#> 5      chr10   chr10_TAD4 1320001 1560000
#> 6      chr10   chr10_TAD5 1560001 1680000

# Gene positions and gene-to-TAD assignment
gene2tadDT_file <- file.path("..", "EXAMPLE", "DATA", "ENCSR489OCU_NCI-H460_all_genes_positions.txt")
gene2tad_dt <- read.delim(gene2tadDT_file,
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
head(gene2tad_dt)
#>    entrezID chromosome  start    end     region
#> 1    347688      chr10  92828  95178 chr10_TAD1
#> 2    439945      chr10 125916 132183 chr10_TAD1
#> 3     10771      chr10 180405 300577 chr10_TAD1
#> 4 100421369      chr10 193770 194736 chr10_TAD1
#> 5     22982      chr10 320130 735621 chr10_TAD2
#> 6 100847086      chr10 687629 687718 chr10_TAD1
```


### Gene expression data

<ul>
<li> RData format expected (i.e. a format that can be imported with <em>load()</em>)</li>
<li>row names are expected to be gene IDs</li>
<li>column names are expected to be sample IDs</li>
<li>can contain extra columns (only those corresponding to sample1 and sample2 IDs will be retained)</li>
</ul>


```r
# used to create classes of gene expression (step 5fc)
rna_fpkmDT_file <- file.path("..", "EXAMPLE", "DATA", "TCGAluad_norm_luad_fpkmDT.RData")
rna_fpkmDT <- get(load(rna_fpkmDT_file))
rna_fpkmDT[1:5,1:5]
#>           TCGA-38-4625-11 TCGA-38-4626-11 TCGA-38-4627-11 TCGA-38-4632-11
#> 100130426    0.000000e+00    0.000000e+00    0.000000e+00    0.000000e+00
#> 100133144    2.109895e-07    7.875051e-08   4.697807e-160    4.436309e-07
#> 100134869    4.146888e-07    3.540477e-07    3.075534e-07    3.535838e-07
#> 10357        1.553667e-05    1.598769e-05    1.860075e-05    1.322551e-05
#> 10431        3.599472e-05    4.673789e-05    5.279937e-05    4.914794e-05
#>           TCGA-44-2655-11
#> 100130426    0.000000e+00
#> 100133144    3.952489e-07
#> 100134869    4.409371e-07
#> 10357        1.804483e-05
#> 10431        3.341402e-05

# used for gene-level differential expression analysis (step)
rnaseqDT_file <- file.path("..", "EXAMPLE", "DATA", "TCGAluad_norm_luad_rnaseqDT_v2.RData")
rnaseqDT <- get(load(rnaseqDT_file))
rnaseqDT[1:5,1:5]
#>           TCGA-38-4625-11 TCGA-38-4626-11 TCGA-38-4627-11 TCGA-38-4632-11
#> 100130426          0.0000          0.0000          0.0000          0.0000
#> 100133144          4.0315          1.4705          0.0000          8.2678
#> 100134869         10.9554          9.1110          7.5432          9.0772
#> 10357            128.3177        129.6157        144.8024        108.5800
#> 10431            659.4230        835.5713        900.8730        878.8159
#>           TCGA-44-2655-11
#> 100130426          0.0000
#> 100133144          5.3467
#> 100134869          8.2126
#> 10357            107.5162
#> 10431            433.5901
```


### Sample information

<ul>
<li>RData format expected (i.e. a format that can be imported with <em>load()</em>)</li>
<li>should be contained in the columns of the gene expression data</li>
<li><em>sample1_file</em> should store sample IDs for <em>cond1</em> and <em>sample2_file</em> those of <em>cond2</em></li>
</ul>


```r
# samples that correspond to cond1 <- <CONDITION_1>
sample1_file <- file.path("..", "EXAMPLE", "DATA", "norm_ID.RData")
samp1 <- get(load(sample1_file))
head(samp1)
#> [1] "TCGA-38-4625-11" "TCGA-38-4626-11" "TCGA-38-4627-11" "TCGA-38-4632-11"
#> [5] "TCGA-44-2655-11" "TCGA-44-2657-11"
# samples that correspond to cond2 <- <CONDITION_2>
sample2_file <- file.path("..", "EXAMPLE", "DATA", "luad_ID.RData")
samp2 <- get(load(sample2_file))
head(samp2)
#> [1] "TCGA-05-4244-01" "TCGA-05-4249-01" "TCGA-05-4250-01" "TCGA-05-4382-01"
#> [5] "TCGA-05-4384-01" "TCGA-05-4389-01"
```



### Import correlation permutation values 

For step 9, case of a provided file:

<ul>
<li>RData format expected (i.e. a format that can be imported with <em>load()</em>)</li>
<li>this is a list of list
<li>each element of the outer list should correspond to permutation value for a given TAD</li>
<li>in the inner list, there should be one element named <em>meanCorr</em></li>
</ul>


```r
all_permutCorr_data <- get(load(file.path("all_sample_corrValues.RData")))
str(all_permutCorr_data[[1]])
#>  num [1:1713] 0.00305 0.06268 0.11925 0.09171 0.36485 ...
```







