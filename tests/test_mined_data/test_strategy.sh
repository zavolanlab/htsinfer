#!/bin/bash

conda activate htsinfer

# Fetching and arranging the data, converting SRA to FASTQ
conda run -n zarp snakemake --snakefile="zarp/workflow/rules/sra_download.smk" \
                            --profile="zarp/profiles/local-conda" \
                            --config samples="mined_test_data_SE.tsv" \
                              outdir="results_sra_downloads" \
                              samples_out="results_sra_downloads/mined_test_data_SE.out.tsv" \
                              log_dir="zarp/logs" \
                              cluster_log_dir="zarp/logs/cluster_log"

# Running HTSinfer over all files

for SAMPLE in results_sra_downloads/*/*.fastq.gz; 
do
htsinfer --output-directory results_htsinfer \
         --cleanup-regime KEEP_ALL \
         --verbosity DEBUG \
         $SAMPLE > results_htsinfer/$(basename "${SAMPLE%%.*}")_result.json
done

# Compilation of the results in a matrix format
python test_strategy.py

# Comparing the results with the mined data, analysis of comparison results
