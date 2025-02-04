import pandas as pd
import numpy as np
from scipy import stats

def single_cell(filepath,min_gene_counts=10,min_sample_counts=10,min_expression=10,min_variance=0.5,p_value_threshold=0.05,genes_per_set=3,output_dir='output/'):
    data = pd.read_csv(filepath)

    data_numeric = data.select_dtypes(include=[np.number])
    data_non_numeric = data.select_dtypes(exclude=[np.number])

    gene_counts = (data_numeric > 0).sum(axis=1)
    genes_to_keep = gene_counts >= min_gene_counts

    sample_counts = data_numeric.sum(axis=0)
    samples_to_keep = sample_counts >= min_sample_counts

    qc_data_numeric = data_numeric.loc[genes_to_keep, samples_to_keep]

    qc_data = pd.concat([data_non_numeric,qc_data_numeric], axis=1) 

    data_numeric = qc_data.select_dtypes(include=[np.number])
    data_non_numeric = qc_data.select_dtypes(exclude=[np.number])

    genes_above_min_expression = (data_numeric >= min_expression).any(axis=1)
    gene_variance = data_numeric.var(axis=1)
    genes_above_min_variance = gene_variance >= min_variance

    filter_condition = genes_above_min_expression & genes_above_min_variance

    filtered_data_numeric = data_numeric[filter_condition].dropna()

    filtered_data_non_numeric = data_non_numeric.loc[filter_condition]

    filtered_data = pd.concat([filtered_data_non_numeric,filtered_data_numeric], axis=1)

    data_numeric = filtered_data.select_dtypes(include=[np.number])
    non_numeric_data = filtered_data.select_dtypes(exclude=[np.number])

    library_sizes = data_numeric.sum(axis=0)
    effective_library_size = np.exp(np.mean(np.log(library_sizes)))
    normalization_factors = library_sizes / effective_library_size

    normalized_data_with_non_numeric = data_numeric.div(normalization_factors, axis=1)

    normalized_data = pd.concat([non_numeric_data,normalized_data_with_non_numeric], axis=1)

    data_numeric = normalized_data.select_dtypes(include=[np.number])
    non_numeric_data = normalized_data.select_dtypes(exclude=[np.number])

    num_columns = len(data_numeric.columns)
    group_A_columns = data_numeric.columns[:num_columns // 2]
    group_B_columns = data_numeric.columns[num_columns // 2:]

    t_test_results = {}
    for gene in data_numeric.index:
        t_statistic, p_value = stats.ttest_ind(data_numeric.loc[gene, group_A_columns],
                                                data_numeric.loc[gene, group_B_columns])
        t_test_results[gene] = {'t_statistic': t_statistic, 'p_value': p_value}

    t_test_df = pd.DataFrame.from_dict(t_test_results, orient='index')
    combined_data = pd.concat([non_numeric_data, t_test_df], axis=1, join='inner')

    # Filter significant genes based on p-value threshold
    significant_genes_filtered = combined_data[combined_data['p_value'] < p_value_threshold]

    # Sort the filtered significant genes
    combined_data_sorted = significant_genes_filtered.sort_values(by='p_value')

    significant_genes_non_numeric = {}
    for gene, result in t_test_results.items():
        p_value = result['p_value']
        if isinstance(p_value, np.ndarray):
            p_value = p_value[0]
        if p_value < p_value_threshold:
            significant_genes_non_numeric[gene] = result

    significant_genes = {}

    for gene, result in significant_genes_non_numeric.items():
        non_numeric_row = non_numeric_data.loc[gene]
    
        result_with_non_numeric = result.copy() 
        result_with_non_numeric.update(non_numeric_row)
    
        significant_genes[gene] = result_with_non_numeric

    sorted_significant_genes = dict(sorted(significant_genes.items(), key=lambda x: x[1]['p_value']))
    sorted_significant_genes_df = pd.DataFrame.from_dict(sorted_significant_genes, orient='index')
    gene_name_column = sorted_significant_genes_df.columns[2]
    significant_gene_names = list(sorted_significant_genes_df[gene_name_column])

    gene_sets = {}
    current_set_index = 1
    for i in range(0, len(significant_gene_names), genes_per_set):
        gene_set_name = f'Gene_Set_{current_set_index}'
        gene_sets[gene_set_name] = significant_gene_names[i:i + genes_per_set]
        current_set_index += 1
        
    gene_name_column=[]
    gene_name_column = combined_data_sorted.columns[0]
    gene_ranks = {}
    for index, row in combined_data_sorted.iterrows():
        gene_name = row[gene_name_column]
        t_statistic = row['t_statistic']
        gene_ranks[gene_name] = t_statistic

    gsea_results = {}
    for gene_set_name, gene_set in gene_sets.items():
        enrichment_scores = [gene_ranks[gene] if gene in gene_ranks else 0 for gene in gene_set]
        enrichment_score = sum(enrichment_scores)
        gsea_results[gene_set_name] = enrichment_score

    sorted_gsea_results = sorted(gsea_results.items(), key=lambda x: np.any(x[1]), reverse=True)
    gsea_results_df = pd.DataFrame(sorted_gsea_results)
    return qc_data,filtered_data,normalized_data,combined_data_sorted,gsea_results_df


