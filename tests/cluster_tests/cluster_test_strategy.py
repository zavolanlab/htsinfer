"""Testing strategy for HTSinfer."""

import json
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from htsinfer.models import Results

ZARP_CMD = 'TEST_PATH=tests/cluster_tests; \
            conda run -n zarp snakemake \
            --snakefile="$TEST_PATH/zarp/workflow/rules/sra_download.smk" \
            --profile="$TEST_PATH/zarp/profiles/local-conda" \
            --config samples="$TEST_PATH/mined_test_data.tsv" \
            outdir="$TEST_PATH/results_sra_downloads" \
            samples_out="$TEST_PATH/results_sra_downloads/mined_data.out.tsv" \
            log_dir="$TEST_PATH/zarp/logs" \
            cluster_log_dir="$TEST_PATH/zarp/logs/cluster_log"'

HTS_CMD = 'TEST_PATH=tests/cluster_tests; \
           conda run -n htsinfer htsinfer \
           --output-directory $TEST_PATH/results_htsinfer \
           --cleanup-regime DEFAULT \
           --verbosity DEBUG \
           --threads 7 \
           --records 200000'

result_data = {}

test_data = Path('tests/cluster_tests/mined_test_data.tsv')
results_sra = Path('tests/cluster_tests/results_sra_downloads/')
results_hts = Path('tests/cluster_tests/results_htsinfer/')

with open(test_data, encoding="utf-8") as tsv_file:
    source = pd.read_csv(tsv_file, sep='\t')
    subprocess.Popen(ZARP_CMD, shell=True,
                     executable='/bin/bash').communicate()

    for index, row in source.iterrows():
        if row['layout'] == 'SE':
            sample_se = results_sra + row['sample'] \
                + '/' + row['sample'] + '.fastq.gz'
            subprocess.Popen(HTS_CMD + ' ' + sample_se
                             + '>' + results_hts +
                             + row['sample'] + '_result.json',
                             shell=True, executable='/bin/bash').communicate()
        else:
            sample_1 = results_sra + row['sample'] \
                + '/' + row['sample'] + '_1.fastq.gz'
            sample_2 = results_sra + row['sample'] \
                + '/' + row['sample'] + '_2.fastq.gz'
            subprocess.Popen(HTS_CMD + ' ' + sample_1 + ' ' + sample_2
                             + '>' + results_hts
                             + row['sample'] + '_result.json',
                             shell=True, executable='/bin/bash').communicate()

    for index, row in source.iterrows():
        if row['layout'] == 'SE':
            with open(results_hts
                      + row['sample'] + '_result.json',
                      encoding="utf-8") as json_file:
                data_model = Results(**json.load(json_file))
                result_data[index] = {
                    'pred_org': data_model.library_source.file_1.short_name,
                    'pred_orient': data_model.read_orientation.file_1,
                    'pred_adapter': data_model.read_layout.file_1.adapt_3,
                    'pred_length_min_1':
                    data_model.library_stats.file_1.read_length.min,
                    'pred_length_max_1':
                    data_model.library_stats.file_1.read_length.max}
                result = source.fillna(
                    pd.DataFrame.from_dict(result_data, orient='index'))
        else:
            with open(results_hts
                      + row['sample']+'_result.json',
                      encoding="utf-8") as json_file:
                data_model = Results(**json.load(json_file))
                result_data[index] = {
                    'pred_org': data_model.library_source.file_1.short_name,
                    'pred_orient': data_model.read_orientation.file_1,
                    'pred_adapter': data_model.read_layout.file_1.adapt_3,
                    'pred_length_min_1':
                    data_model.library_stats.file_1.read_length.min,
                    'pred_length_max_1':
                    data_model.library_stats.file_1.read_length.max,
                    'pred_length_min_2':
                    data_model.library_stats.file_2.read_length.min,
                    'pred_length_max_2':
                    data_model.library_stats.file_2.read_length.max}
                result = source.fillna(
                    pd.DataFrame.from_dict(result_data, orient='index'))

    # Comparison of results
    result['match_org'] = np.where(
        result['org'] == result['pred_org'], True, False)
    result['match_length'] = np.where(
        result['pred_length_max_1'] == result['length_max_1'], True, False)
    result['match_adapter'] = result.apply(
        lambda x: str(x.pred_adapter) in str(x.adapter), axis=1)
    print(result)
    pd.DataFrame.to_csv(
        result, 'tests/cluster_tests/mined_test_data_result.tsv',
        sep='\t', index=False)
