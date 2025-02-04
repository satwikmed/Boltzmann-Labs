import os
import boto3

def download_data(qc_data, filtered_data, normalized_data, combined_data_sorted, gsea_results_df, output_directory,experiment_id,job_name):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    qc_data_file = os.path.join(output_directory, 'qc_data.csv')
    filtered_data_file = os.path.join(output_directory, 'filtered_data.csv')
    normalized_data_file = os.path.join(output_directory, 'normalized_data.csv')
    combined_data_sorted_file = os.path.join(output_directory, 'combined_data_sorted.csv')
    gsea_results_file = os.path.join(output_directory, 'gsea_results.csv')

    qc_data.to_csv(qc_data_file, index=False)
    filtered_data.to_csv(filtered_data_file, index=False)
    normalized_data.to_csv(normalized_data_file, index=False)
    combined_data_sorted.to_csv(combined_data_sorted_file, index=False)
    gsea_results_df.to_csv(gsea_results_file, index=False)

    print(f"Data downloaded successfully to directory: {output_directory}/{experiment_id}/{job_name}")

def upload_to_s3(output_directory, bucket_name, object_key,experiment_id,job_name):
    try:
        s3_client = boto3.client('s3')

        for root, dirs, files in os.walk(output_directory):
            for file in files:
                local_file_path = os.path.join(root, file)

                object_key_file = os.path.join(object_key, file).replace("\\", "/")
                s3_client.upload_file(local_file_path, bucket_name, object_key_file)
                print(f"File uploaded successfully to S3 bucket: {bucket_name}/{object_key_file}")
        return True
    except Exception as e:
        print(f"Error uploading files to S3: {e}")
        return False
    
def download_file(bucket_name, object_name, local_file_path):
    s3 = boto3.client('s3')
    try:
        s3.download_file(bucket_name, object_name, local_file_path)
        print(f"File downloaded from {bucket_name}/{object_name} to {local_file_path}")
        return local_file_path
    except Exception as e:
        print(f"Error downloading file: {e}")
