"""Testing strategy for HTSinfer."""

import json
import subprocess as sp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from htsinfer.models import Results

records = 100000
min_match = 2
min_freq = 0.3

ZARP_CMD = 'TEST_PATH=tests/cluster_tests; \
            ZARP_PATH=../zarp; \
            conda run -n zarp snakemake --unlock \
            --snakefile="$ZARP_PATH/workflow/rules/sra_download.smk" \
            --profile="$ZARP_PATH/profiles/local-conda" \
            --config samples="$TEST_PATH/mined_test_data.tsv" \
            outdir="$TEST_PATH/results_sra_downloads" \
            samples_out="$TEST_PATH/results_sra_downloads/mined_data.out.tsv" \
            log_dir="$ZARP_PATH/logs" \
            cluster_log_dir="$ZARP_PATH/logs/cluster_log"'

FOLDER_CMD = 'FOLDER=$(date +%m%d_%H%M%S); \
              TEST_PATH=tests/cluster_tests; \
              mkdir $TEST_PATH/results_htsinfer/$FOLDER; \
              echo "$FOLDER" > $TEST_PATH/results_htsinfer/temp.txt'

HTS_CMD = 'TEST_PATH=tests/cluster_tests; \
           FOLDER=`cat $TEST_PATH/results_htsinfer/temp.txt`; \
           conda run -n htsinfer htsinfer \
           --output-directory $TEST_PATH/results_htsinfer/$FOLDER \
           --cleanup-regime KEEP_ALL \
           --verbosity DEBUG \
           --read-layout-min-match-percentage ' + str(min_match) + '\
           --read-layout-min-frequency-ratio ' + str(min_freq) + '\
           --threads 7 \
           --records ' + str(records)

result_data = {}
graph_data = {}

# test files
PACKAGE_DIR = Path(__file__).resolve().parents[3] / "htsinfer"
ZARP_DIR = Path(__file__).resolve().parent / "zarp"
RESULTS_SRA_DIR = Path(__file__).resolve().parent / "results_sra_downloads"
RESULTS_HTS_DIR = Path(__file__).resolve().parent / "results_htsinfer"
MINED_DATA = Path(__file__).resolve().parent / "mined_test_data.tsv"


class TestDownloadSRASamples:
    """Download the SRA Samples."""

    def test_zarp_download(self):
        sp.Popen(ZARP_CMD, shell=True,
                 executable='/bin/bash').communicate()


class TestHTSinfer:
    """Test HTSinfer on downloaded samples."""

    def test_htsinfer_se_pe(self):
        with open(MINED_DATA, encoding="utf-8") as tsv_file:
            source = pd.read_csv(tsv_file, sep='\t')
            sp.Popen(FOLDER_CMD, shell=True,
                     executable='/bin/bash').communicate()
            for index, row in source.iterrows():
                if row['layout'] == 'SE':
                    sample_se = str(RESULTS_SRA_DIR) + '/' + row['sample'] \
                                + '/' + row['sample'] + '.fastq.gz'
                    sp.Popen(HTS_CMD + ' ' + sample_se + ' > '
                             + str(RESULTS_HTS_DIR) + '/$FOLDER/'
                             + row['sample'] + '_result.json',
                             shell=True,
                             executable='/bin/bash').communicate()
                else:
                    sample_1 = str(RESULTS_SRA_DIR) + '/' + row['sample'] \
                        + '/' + row['sample'] + '_1.fastq.gz'
                    sample_2 = str(RESULTS_SRA_DIR) + '/' + row['sample'] \
                        + '/' + row['sample'] + '_2.fastq.gz'
                    sp.Popen(HTS_CMD + ' ' + sample_1 + ' '
                             + sample_2 + ' > ' + str(RESULTS_HTS_DIR)
                             + '/$FOLDER/' + row['sample'] + '_result.json',
                             shell=True,
                             executable='/bin/bash').communicate()

    def test_write_results(self):
        with open(MINED_DATA, encoding="utf-8") as tsv_file:
            source = pd.read_csv(tsv_file, sep='\t')
            results_folder = (
                sp.check_output('TEST_P=tests/cluster_tests/results_htsinfer; \
                                FOLDER=`cat $TEST_P/temp.txt`; \
                                echo $FOLDER', shell=True,
                                executable='/bin/bash').strip().decode(
                                'ascii'))
            for index, row in source.iterrows():
                if row['layout'] == 'SE':
                    with open(str(RESULTS_HTS_DIR) + '/' + results_folder + '/'
                              + row['sample'] + '_result.json',
                              encoding="utf-8") as json_file:
                        data_model = Results(**json.load(json_file))
                        result_data[index] = {
                            'pred_org':
                            data_model.library_source.file_1.short_name,
                            'pred_orient':
                            data_model.read_orientation.file_1,
                            'pred_adapter':
                            data_model.read_layout.file_1.adapt_3,
                            'pred_length_min_1':
                            data_model.library_stats.file_1.read_length.min,
                            'pred_length_max_1':
                            data_model.library_stats.file_1.read_length.max}
                        result = source.fillna(
                            pd.DataFrame.from_dict(result_data,
                                                   orient='index'))
                else:
                    with open(str(RESULTS_HTS_DIR) + '/' + results_folder + '/'
                              + row['sample'] + '_result.json',
                              encoding="utf-8") as json_file:
                        data_model = Results(**json.load(json_file))
                        result_data[index] = {
                            'pred_org':
                            data_model.library_source.file_1.short_name,
                            'pred_orient':
                            data_model.read_orientation.file_1,
                            'pred_adapter':
                            data_model.read_layout.file_1.adapt_3,
                            'pred_length_min_1':
                            data_model.library_stats.file_1.read_length.min,
                            'pred_length_max_1':
                            data_model.library_stats.file_1.read_length.max,
                            'pred_length_min_2':
                            data_model.library_stats.file_2.read_length.min,
                            'pred_length_max_2':
                            data_model.library_stats.file_2.read_length.max}
                        result = source.fillna(
                            pd.DataFrame.from_dict(result_data,
                                                   orient='index'))

            # Comparison of results
            result['match_org'] = np.where(
                result['org'] == result['pred_org'], True, False)
            result['match_length'] = np.where(
                result['pred_length_max_1'] ==
                result['length_max_1'], True, False)
            result['match_adapter'] = result.apply(
                lambda x: str(x.pred_adapter) in str(x.adapter), axis=1)
            pd.DataFrame.to_csv(
                result, 'tests/cluster_tests/results_htsinfer/'
                + results_folder + '/mined_test_data_result.tsv',
                sep='\t', index=False)

        with open(MINED_DATA, encoding="utf-8") as tsv_file:
            source = pd.read_csv(tsv_file, sep='\t')
            results_folder = (
                sp.check_output('TEST_P=tests/cluster_tests/results_htsinfer; \
                                FOLDER=`cat $TEST_P/temp.txt`; \
                                echo $FOLDER', shell=True,
                                executable='/bin/bash').strip().decode(
                                'ascii'))
            for index, row in source.iterrows():
                if row['layout'] == 'SE':
                    with open(str(RESULTS_HTS_DIR) + '/' + results_folder + '/'
                              + 'read_layout_' + row['sample'] + '.fastq.json',
                              encoding="utf-8") as json_file:
                        graph_model = json.load(json_file)
                        graph_data[index] = {
                            '1_adapt_1': graph_model['data'][0][0],
                            '1_percent_1': graph_model['data'][0][1],
                            '1_adapt_2': graph_model['data'][1][0],
                            '1_percent_2': graph_model['data'][1][1]}
                else:
                    with open(str(RESULTS_HTS_DIR) + '/' + results_folder + '/'
                              + 'read_layout_' + row['sample']
                              + '_1.fastq.json',
                              encoding="utf-8") as json_file_1, open(
                              str(RESULTS_HTS_DIR) + '/' + results_folder + '/'
                              + 'read_layout_' + row['sample'] +
                              '_2.fastq.json',
                              encoding="utf-8") as json_file_2:
                        graph_model_1 = json.load(json_file_1)
                        graph_model_2 = json.load(json_file_2)
                        graph_data[index] = {
                            '1_adapt_1': graph_model_1['data'][0][0],
                            '1_percent_1': graph_model_1['data'][0][1],
                            '1_adapt_2': graph_model_1['data'][1][0],
                            '1_percent_2': graph_model_1['data'][1][1],
                            '2_adapt_1': graph_model_2['data'][0][0],
                            '2_percent_1': graph_model_2['data'][0][1],
                            '2_adapt_2': graph_model_2['data'][1][0],
                            '2_percent_2': graph_model_2['data'][1][1]}
            graph_result = pd.DataFrame.from_dict(graph_data, orient='index')
            result_final = pd.concat([result, graph_result],
                                     axis=1, join="inner")
            # print(result_final)
            pd.DataFrame.to_csv(result_final,
                                'tests/cluster_tests/results_htsinfer/'
                                + results_folder + '/graph_result.tsv',
                                sep='\t', index=False)

            # Barplot
            result_final.loc[result_final['pred_adapter'].notnull(),
                             'pred_adapter'] = 1
            result_final['pred_adapter'] = (
                result_final['pred_adapter'].fillna(0))
            fig, axs = plt.subplots(1, figsize=[8, 8])
            ax = sns.barplot(x='pred_adapter', y='pred_adapter',
                             data=result_final,
                             estimator=lambda x: len(x) /
                             len(result_final) * 100)
            ax.set(ylabel="Percent")
            ax.bar_label(ax.containers[0], fmt='%.f%%')
            ax.set(xticklabels=["Not identified", "Identified"])
            ax.set(title='Number and fraction of identified adapters\nID: '
                   + str(results_folder)
                   + '\nNo. of records: ' + str(records)
                   + '\nRead layout min-match percentage: ' + str(min_match)
                   + '\nRead layout min frequency ratio: ' + str(min_freq))
            plt.savefig(str(RESULTS_HTS_DIR) + '/' + str(results_folder)
                        + '/1_Barplot_predicted_adapters', dpi=100)

            # Histogram of 1st predicted adapter percent #1
            # # Drop 0 percentages of both SE and PE reads
            result_final_n = result_final.drop(result_final[
                (result_final['1_percent_1'] == 0) &
                (result_final['2_percent_1'] == 0) |
                (result_final['1_percent_1'] == 0) &
                result_final['2_percent_1'].isna()].index)
            all_percent = pd.concat([result_final_n['1_percent_1'],
                                     result_final_n['2_percent_1']])
            all_percent = all_percent[all_percent != 0]
            fig, axs = plt.subplots(1, figsize=[8, 8])
            sns.histplot(data=all_percent, binwidth=2).set(
                title='Fraction of reads containing most '
                      + 'prevalent adapter\nID: ' + str(results_folder)
                      + '\nNo. of records: ' + str(records)
                      + '\nRead layout min-match percentage: ' + str(min_match)
                      + '\nRead layout min frequency ratio: ' + str(min_freq))
            plt.xlim(0, 100)
            plt.savefig(str(RESULTS_HTS_DIR) + '/' + str(results_folder)
                        + '/2_Hist_1st_pred_adapter_full.png', dpi=100)

            # Histogram of 1st predicted adapter percent #2
            fig, axs = plt.subplots(1, figsize=[8, 8])
            sns.histplot(data=all_percent, binwidth=0.2).set(
                title='Fraction of reads containing most '
                      + 'prevalent adapter\nID: '
                      + str(results_folder)
                      + '\nNo. of records: ' + str(records)
                      + '\nRead layout min-match percentage: ' + str(min_match)
                      + '\nRead layout min frequency ratio: ' + str(min_freq))
            plt.xlim(0, 10)
            plt.savefig(str(RESULTS_HTS_DIR) + '/' + str(results_folder)
                        + '/3_Hist_1st_pred_adapter_10.png', dpi=100)

            # Histogram of 1st vs 2nd predicted adapter ratio
            result_final_n = result_final.drop(result_final[
                (result_final['1_percent_1'] == 0) &
                (result_final['2_percent_1'] == 0) |
                (result_final['1_percent_1'] == 0) &
                result_final['2_percent_1'].isna()].index)
            result_final_n['1_ratio'] = (
                result_final_n['1_percent_1'] / (
                    result_final_n['1_percent_2'] + 0.01))
            result_final_n['2_ratio'] = (
                result_final_n['2_percent_1'] / (
                    result_final_n['2_percent_2'] + 0.01))
            all_ratios = pd.concat(
                [result_final_n['1_ratio'], result_final_n['2_ratio']])
            fig, axs = plt.subplots(1, figsize=[8, 8])
            sns.histplot(data=all_ratios).set(
                title='Fraction of reads with most prevalent adapter '
                + 'vs. second most prevalent\nID: '
                + str(results_folder)
                + '\nNo. of records: ' + str(records)
                + '\nRead layout min-match percentage: ' + str(min_match)
                + '\nRead layout min frequency ratio: ' + str(min_freq),
                xscale="log")
            plt.savefig(str(RESULTS_HTS_DIR) + '/' + str(results_folder)
                        + '/4_Hist_1st_vs_2nd_pred_adapter_ratio.png', dpi=100)
