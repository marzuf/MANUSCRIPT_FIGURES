##### Main settings

Set main settings in "main_settings.R" (cf. template file: "TEMPLATE_setting_file.R").

Do not change the name of this file.

The current "main_settings.R" file was used to run the example as explained below.


##### Additional settings

Set additional (e.g. dataset-specific) settings in another file (cf. template file: "example_setting_file.R")

The name of this file can be user-defined.

Settings from this file will overwrite settings from "main_settings.R" (hence can be used for this purpose).


##### Launch the pipeline

```{bash}
 ./run_pipeline.sh <setting_file> <stepNbr>
```

`setting_file` should be path to the additional settings file.

`stepNbr` indicates the step of the pipeline that should be run (should be what comes before the "_" in the script name, e.g. 10 or 10provided for the 10-th step)

The names of the scripts explicitly describe to which step they correspond, namely:
<ul>
<li><em>1_prepGeneData.R</em>: prepare some data for the rest of the pipeline</li>
<li>*2_runGeneDE.R*: gene-level differential expression analysis (needed for logFC values)</li>
<li>*3_runMeanTADLogFC.R*: prepare TAD-level average FC</li>
<li>*4_runMeanTADCorr.R*: prepare TAD-level average correlation</li>
<li>*5fc_runPermutationsMedian.R*: run permutation for logFC values (can take some times depending on the number of permutations !)</li>
<li>*5corr_runPermutationsCorr.R*: run permutation for correlation values</li>
<li>*6_runPermutationsMeanLogFC.R*: compute TAD-level average FC for the permutation data</li>
<li>*7_runPermutationsMeanTADCorr.R*:  compute TAD-level average correlation for the permutation data</li>
<li>*8_runEmpPvalMeanTADLogFC.R*: compute empirical p-values for the FC</li>
<li>*9_runEmpPvalMeanTADCorr.R*: compute empirical p-values for correlation</li>
<li>*10_runEmpPvalCombined.R*: combine empirical p-values for FC and correlation</li>
</ul>

##### Outputs

Each script outputs files in a folder with the same name as the script (e.g. files written in "1_prepGeneData" for "1_prepGeneData.R").


##### Dependencies

[pigz](https://zlib.net/pigz) is required for fast saving in Rdata files in the scripts 5fc and 6 (we used version 2.4).


R packages used in the pipeline: foreach, doMC, dplyr


##### Notes

The pipeline is designed to run the steps one after the other (some depedencies exist among the scripts, e.g. the 11th requires the 10th and 9th to have been run).

**WARNING**: 
In step 9, there are two possibilities for retrieving correlation permutation data to be set in the setting files:

<ul>
<li>if `all_permutCorr_data` is a folder:</li>
<ul>
<li>the correlation values from permutation data are loaded from files in `all_permutCorr_data` that (recursively) match the pattern "corrMatchPattern" (set to "meanCorr_sample_around_TADs_sameNbr.Rdata" in the main settings but can be overwritten in the additional settings)</li>
<li>if "refineMatchPattern" is provided, used for a second pattern matching to refine file retrieval</li>
<li>if "corrDiscardPattern" is provided in setting files, files that match that pattern are discarded</li>
<li>if "nbrCorrPermutCheck" is provided, it is checked that "nbrCorrPermutCheck" files have been loaded</li>
<li>it is expected to correspond to the format of script 7: each file should correspond to a list of lists storing a correlation value in my_list[[idx]][["meanCorr"]] with as many "idx" as TADs
</ul>
<li>if `all_permutCorr_data` is a file:</li>
<ul>
<li>`all_permutCorr_file` should provide the path to correlation values for the permutation data (we provide the values from our permuation data in "data/all_sample_corrValues.Rdata"); the data should be a Rdata file containing a list/vector of permutation correlation values</li>
</ul>
</ul>



##### Example

Example for launching the pipeline (using provided data for intra-TAD correlation):

```{bash}
 ./example_run_pipeline.sh example_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAlusc_norm_lusc.R 1 2 3 4 5fc 6 8 9provided 10provided
```

Example for launching the pipeline (full pipeline)

```{bash}
 ./example_run_pipeline.sh example_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAlusc_norm_lusc.R 1 2 3 4 5corr 5fc 6 7 8 9 10
```
