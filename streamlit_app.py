import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import boto3

# Read AWS credentials from environment variables
@st.cache_data
def get_adata():
    aws_access_key_id = st.secrets['AWS_ACCESS_KEY_ID']
    aws_secret_access_key = st.secrets['AWS_SECRET_ACCESS_KEY']

    # Initialize S3 client
    s3 = boto3.client('s3', aws_access_key_id=aws_access_key_id, aws_secret_access_key=aws_secret_access_key)

    bucket_name = 'testzl57208'
    file_key = 'poplar_axl_bud_A_adata.h5ad'
    local_file_name = 'poplar_axl_bud_A_adata.h5ad'
    # print(f'Reading file: {local_file_name}')

    # Download file from S3
    s3.download_file(bucket_name, file_key, local_file_name)

    # Read the AnnData file using Scanpy
    adata = sc.read(local_file_name)
    return adata

## main page
st.markdown("# Data Evaluation App")

st.write("We are so glad to see you here. ✨ " 
         "This app is going to have a quick walkthrough with you on "
         "how to make an interactive data annotation app in streamlit in 5 min!")

def main():
    ## prepare the data
    # gene_ids = adata.var.index.tolist()
    gene_ids = ['test','ann1.Glyma.02G228100', 'ann1.Glyma.15G127900']
    
    ## Define sidebar #####################################
    st.sidebar.header('Plot Configuration')
    st.write('Please select the plot type (UMAP or spatial), and gene ID for the plot')
    # plot_type = st.sidebar.selectbox('Select plot type', ['UMAP', 'Spatial'])
    plot_type = st.sidebar.radio('Select plot type', ['UMAP', 'Spatial'], horizontal=True)
    gene_name = st.sidebar.selectbox('Enter gene name for expression plot', ['', *gene_ids])
    
    
    
    
    
    st.write("Here we are at the end of getting started with streamlit! Happy Streamlit-ing! :balloon:")

    
    
    

if __name__ == '__main__':
    main()
