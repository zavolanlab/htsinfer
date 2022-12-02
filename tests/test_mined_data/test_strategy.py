import json
import pandas as pd

result_data = {}

with open('mined_test_data_SE.tsv') as tsv_file:
    source = pd.read_csv(tsv_file, sep='\t')
    for index, row in source.iterrows():
        with open('results_htsinfer/'+row['sample']+'_result.json') as json_file:
            data = json.load(json_file)
            tmp = { 'Pred. Organism':data["library_source"]["file_1"]["short_name"], 
                    'Pred. Read orientation':data["read_orientation"]["file_1"], 
                    'Pred. 3\'-adapter':data["read_layout"]["file_1"]["adapt_3"], 
                    'Pred. Read length':data["library_stats"]["file_1"]["read_length"]}
            result_data[index] = tmp
    source.fillna(pd.DataFrame.from_dict(result_data, orient='index')).to_csv('mined_test_data_SE_result.tsv', sep='\t', index=False)
