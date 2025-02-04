import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import boto3
from sklearn.cluster import KMeans

def save_plot_as_png(plt, filename):
    plt.savefig(filename, format='png')
    plt.close()

def upload_file_to_s3(local_file_path, bucket_name, object_key):
    s3 = boto3.client('s3')

    try:
        s3.upload_file(local_file_path, bucket_name, object_key)
        print(f"File uploaded successfully to S3 bucket: {bucket_name}/{object_key}")
        return True
    except Exception as e:
        print(f"Error uploading file to S3: {e}")
        return False

def visualize_data(filtered_data, combined_data_sorted, gsea_results_df, output_dir, experiment_id, user_id):
    def heatfilter(filtered_data):
        filtered_data = filtered_data.set_index(filtered_data.columns[0])

        plt.figure(figsize=(10, 8))
        sns.heatmap(filtered_data, cmap='coolwarm', robust=True, cbar_kws={'label': 'Expression Level'}, alpha=0.7) # Adjust alpha for transparency
        plt.title('Heatmap of Filtered Data')
        plt.xlabel('Samples')
        plt.ylabel('Genes')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        heatmap_filename = os.path.join(output_dir, 'heatmap.png')
        save_plot_as_png(plt, heatmap_filename)
        upload_file_to_s3(heatmap_filename, bucket_name='boltbio', object_key=f'{user_id}/SingleCell/{experiment_id}/Visualizations/heatmap.png')

    def p_value(sorted_significant_genes_df):
        gene_name_column = combined_data_sorted.columns[0]
        sorted_significant_genes_df[gene_name_column] = sorted_significant_genes_df[gene_name_column].astype(str)

        plt.figure(figsize=(12, 8))
        plt.barh(sorted_significant_genes_df[gene_name_column], sorted_significant_genes_df['p_value'], color='skyblue', alpha=0.7) # Adjust alpha for transparency
        plt.xlabel('Gene Names')
        plt.ylabel('P-Values')
        plt.title('P-values of Significant Genes')
        plt.gca().invert_yaxis()
        plt.grid(axis='x')
        plt.tight_layout()

        p_values_filename = os.path.join(output_dir, 'p_values.png')
        save_plot_as_png(plt, p_values_filename)
        upload_file_to_s3(p_values_filename, bucket_name='boltbio', object_key=f'{user_id}/SingleCell/{experiment_id}/Visualizations/p_values.png')

    def gsea(gsea_results_df):
        gene_sets = [gene_set for gene_set in gsea_results_df[0]]
        enrichment_scores = [score for score in gsea_results_df[1]]

        plt.figure(figsize=(12, 8))
        plt.bar(gene_sets, enrichment_scores, color='skyblue', alpha=0.7) # Adjust alpha for transparency
        plt.xlabel('Gene Sets')
        plt.ylabel('Enrichment Score')
        plt.title('Enrichment Scores of Gene Sets')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        gsea_filename = os.path.join(output_dir, 'gsea.png')
        save_plot_as_png(plt, gsea_filename)
        upload_file_to_s3(gsea_filename, bucket_name='boltbio', object_key=f'{user_id}/SingleCell/{experiment_id}/Visualizations/gsea.png')

    def perform_clustering(filtered_data, num_clusters=3):
        data_numeric = filtered_data.select_dtypes(include=[np.number])

        kmeans = KMeans(n_clusters=num_clusters)
    
        kmeans.fit(data_numeric)
    
        labels = kmeans.labels_
    
        centers = kmeans.cluster_centers_
    
        plt.figure(figsize=(8, 6))
    
        for i in range(num_clusters):
            cluster_points = data_numeric[labels == i]
        
            plt.scatter(cluster_points.iloc[:, 0], cluster_points.iloc[:, 1], label=f'Cluster {i+1}', alpha=0.7) # Adjust alpha for transparency
    
        plt.scatter(centers[:, 0], centers[:, 1], color='black', marker='x', s=100, label='Centroids')
    
        plt.title('Clusters')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
    
        plt.legend()
    
        plt.grid(True)
        plt.tight_layout()

        clustering_filename = os.path.join(output_dir, 'clustering.png')
        save_plot_as_png(plt, clustering_filename)
        upload_file_to_s3(clustering_filename, bucket_name='boltbio', object_key=f'{user_id}/SingleCell/{experiment_id}/Visualizations/clustering.png')


    heatfilter(filtered_data)
    p_value(combined_data_sorted)
    gsea(gsea_results_df)
    perform_clustering(filtered_data)

        
