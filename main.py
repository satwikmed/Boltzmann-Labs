import os
import json
import sys,traceback

bucket_name='boltbio'

if __name__ == '__main__':
    task_name = os.getenv('TASK_NAME')
    user_input = os.getenv('USER_INPUT')

    user_input = json.loads(user_input)

    job_name = user_input['job_name']
    experiment_id = user_input['experiment_id']
    user_id = user_input['user_id']
    csv_aws_s3 = user_input['datapath']
    property_name = user_input['property_name']
    choice=user_input['choice']


    
    if task_name == "Single Cell":

        job_name = user_input['job_name']
        experiment_id = user_input['experiment_id']
        user_id = user_input['user_id']
        csv_aws_s3 = user_input['datapath']
        property_name = user_input['property_name']
        choice=user_input['choice']

        if property_name == "Workflow":
            try:
                from single_cell import single_cell
                from download_upload import download_data,upload_to_s3,download_file
                job_name = user_input['job_name']
                experiment_id = user_input['experiment_id']
                user_id = user_input['user_id']
                csv_aws_s3 = user_input['datapath']
                try:
                    if choice=='true':
                        local_path = f"output/{experiment_id}_{job_name}.csv"
                        filepath=download_file(bucket_name,object_name=csv_aws_s3,local_file_path=local_path)
                        print(filepath)
                        print(f"File downloaded from S3 bucket '{bucket_name}' to '{local_path}'")
                    else:
                        filepath='/Users/satwikmedipalli/single_cell/GSE237873.csv'
                    qc_data,filtered_data,normalized_data,combined_data_sorted,gsea_results_df=single_cell(filepath,min_gene_counts=10,min_sample_counts=10,min_expression=10,min_variance=0.5,p_value_threshold=0.05,genes_per_set=3)
                    download_data(qc_data, filtered_data, normalized_data, combined_data_sorted, gsea_results_df, output_directory='output/',experiment_id=experiment_id,job_name=job_name)
                    upload_to_s3(output_directory='output/',bucket_name='boltbio',object_key=f'{user_id}/Bulk_RNA/{experiment_id}/{property_name}',experiment_id=experiment_id,job_name=job_name)
                except :
                    print(f"Error downloading file from S3")
                    exc = sys.exception()
                    traceback.print_exception(exc, limit=2, file=sys.stdout)
            except KeyboardInterrupt:
                print('Interrupted')
                try:
                    sys.exit(0)
                except SystemExit:
                    os._exit(0)

        
        
        if property_name == "Visualizations":
            from visualizations import visualize_data
            from download_upload import download_file
            try:
                job_name = user_input['job_name']
                experiment_id = user_input['experiment_id']
                user_id = user_input['user_id']
                csv_aws_s3 = user_input['datapath']
                property_name = user_input['property_name']

                try:
                    from single_cell import single_cell
                    from visualizations import visualize_data
                    if choice=='true':
                        local_path = f"visualizations/{experiment_id}_{job_name}.csv"
                        filepath=download_file(bucket_name,object_name=csv_aws_s3,local_file_path=local_path)
                        print(f"File downloaded from S3 bucket '{bucket_name}' to '{local_path}'")
                    else:
                        filepath='/Users/satwikmedipalli/single_cell/GSE237873.csv'
                        output_dir='visualizations/'
                    qc_data,filtered_data,normalized_data,combined_data_sorted,gsea_results_df=single_cell(filepath,min_gene_counts=10,min_sample_counts=10,min_expression=10,min_variance=0.5,p_value_threshold=0.05,genes_per_set=3)
                    visualize_data(filtered_data, combined_data_sorted, gsea_results_df,output_dir='visualizations/',experiment_id=experiment_id,user_id=user_id)
                except :
                    print(f"Error downloading file from S3")
                    exc = sys.exception()
                    traceback.print_exception(exc, limit=2, file=sys.stdout)
            except KeyboardInterrupt:
                print('Interrupted')
                try:
                    sys.exit(0)
                except SystemExit:
                    os._exit(0)
        