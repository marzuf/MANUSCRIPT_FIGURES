#!/usr/bin/bash


#  ./example_run_pipeline.sh example_ENCSR489OCU_NCI-H460_40kb_run_settings_TCGAlusc_norm_lusc.R 11 


scriptFolder="."  # where the script files are located

set -e

args=( "$@" )
args_len=${#args[@]}

settingF=${args[0]}

i=1
while [[ $i -lt args_len ]]; do
    scriptF="$(ls $scriptFolder | grep -P ^${args[$i]}_.+R$)"
    cmd="Rscript $scriptFolder/$scriptF $settingF"
    echo "> $cmd"
    $cmd
    ((i++))
done


