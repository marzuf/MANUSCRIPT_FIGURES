##### Main settings

Set main settings in "main_settings.R"

Do not change the name of this file.

##### Additional settings

Set additional (e.g. dataset-specific) settings in another file (e.g. "example_setting_file.R")

The name of this file can be user-defined.

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
