##### Main settings

Set main settings in "main_settings.R" (cf. template file: "TEMPLATE_setting_file.R").

Do not change the name of this file.


##### Additional settings

Set additional (e.g. dataset-specific) settings in another file (cf. template file: "example_setting_file.R")

The name of this file can be user-defined.

Settings from this file will overwrite settings from "main_settings.R" (hence can be used for this purpose).


##### Launch the pipeline

```{bash}
 ./run_pipeline.sh <setting_file> <stepNbr>
```

`setting_file` should be path to the additional file settings.


`stepNbr` indicates the step of the pipeline that should be run (should be what comes before the "_" in the script name)


##### Outputs

Each script outputs files in a folder with the same name as the script (e.g. files written in "1_prepGeneData" for "1_prepGeneData.R").


##### Dependencies

[pigz](https://zlib.net/pigz) is required for fast saving in Rdata files in the scripts 5fc and 6 (we used version 2.4).


R packages used in the pipeline: foreach, doMC, dplyr


##### Notes

The pipeline is designed to run the steps one after the other (some depedencies exist among the scripts, e.g. the 11th requires the 10th and 9th to have been run).

**WARNING**: 


In step 9 (), the correlation values from permutation are loaded from files that (recursively) match the pattern "meanCorr_sample_around_TADs_sameNbr.Rdata" located in the `dirname(dirname(pipOutFold))` folder (this might need to be changed !) ("pipOutFold" is set in setting file), and files that contain "RANDOM" or "PERMUT" are discarded. 

We also provide correlation values from our permutation data ("data/all_sample_corrValues.Rdata"). To run the pipeline with these provided data, use the scripts 9provided and 10provided (instead of 9 and 10). To use your own provided permutation correlation values, replace the line `all_sampleCorr_files <- "data/all_sample_corrValues.Rdata"` with the path to your data ( `all_sampleCorr_files <- "path_to_data"`; the data should be a Rdata file containing a list/vector of permutation correlation values).


##### Example

Example for launching the pipeline (using provided data for intra-TAD correlation):

```{bash}
 ./example_run_pipeline.sh example_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAlusc_norm_lusc.R 1 2 3 4 5fc 6 8 9provided 10provided
```

Example for launching the pipeline (full pipeline)

```{bash}
 ./example_run_pipeline.sh example_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAlusc_norm_lusc.R 1 2 3 4 5corr 5fc 6 7 8 9 10
```
