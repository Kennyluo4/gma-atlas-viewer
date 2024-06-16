import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import boto3        # for s3 file transfer
# import paramiko     # for sftp file transfer
# from io import StringIO
import os

@st.cache_data
def get_adata_aws(file_name):
    '''read files stored on AWS S3'''
    # Read AWS credentials from environment variables
    aws_access_key_id = st.secrets['AWS_ACCESS_KEY_ID']
    aws_secret_access_key = st.secrets['AWS_SECRET_ACCESS_KEY']

    # Initialize S3 client
    s3 = boto3.client('s3', aws_access_key_id=aws_access_key_id, aws_secret_access_key=aws_secret_access_key)

    bucket_name = 'soybeanatlas'
    file_key = file_name
    local_file_name = file_name
    # print(f'Reading file: {local_file_name}')

    # Download file from S3
    s3.download_file(bucket_name, file_key, local_file_name)

    # Read the AnnData file using Scanpy
    adata = sc.read(local_file_name)
    return adata

def get_adata_sftp(file_name):
    '''read files stored on cluster'''
    hostname = 'xfer.gacrc.uga.edu'
    # st.secrets['xferhost']
    port = 22
    username = ''
    # st.secrets['userid']
    password = ''
    # st.secrets['xferpass']
    remote_path = '/' + file_name
    with open(os.path.expanduser('~/.ssh/id_rsa'), 'r') as file:
        key_content = file.read()
    try:
        private_key = paramiko.RSAKey.from_private_key(StringIO(key_content), password='')

        transport = paramiko.Transport((hostname, port))
        transport.connect(username=username, pkey=private_key)

        sftp = paramiko.SFTPClient.from_transport(transport)
        with sftp.file(remote_path, 'r') as remote_file:
            file_content = remote_file.read()
        sftp.close()
        transport.close()
        return file_content
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return None

def main():
    ## dic for matching input sample same and stored adata file name
    sample_dic = {'heart stage seed (spatial)': 'gma_sp_HSA_fromSeurat.h5ad','cotyledon stage seed (spatial)':'gma_sp_CS2A_fromSeurat.h5ad', 'early maturation stage seed (spatial)':'gma_sp_ESA_fromSeurat.h5ad',
                'heart stage seed (snRNA)': 'gma_sp_CS2A_fromSeurat.h5ad', 'cotyledon stage seed (snRNA)':'sc_adt_CS.h5ad'}
    
    ###############################################
    ##                 Sidebar                  ###
    ###############################################
    # with st.sidebar.form():           # use if want to add a Submit button before change everything
    
    st.sidebar.header('Plot Configuration')
    st.sidebar.markdown('## Please select a dataset:')
    lib_type = st.sidebar.selectbox('Library type', [ '---Please choose---','snRNA-seq', 'spRNA-seq'])
    
    if lib_type == 'snRNA-seq':
        sample_name = st.sidebar.selectbox('Tissue', ['---Please choose---', 'heart stage seed (snRNA)', 'cotyledon stage seed (snRNA)'])
    elif lib_type == 'spRNA-seq':
        sample_name = st.sidebar.selectbox('Tissue', ['---Please choose---', 'heart stage seed (spatial)', 'cotyledon stage seed (spatial)', 'early maturation stage seed (spatial)'])
    else:
        sample_name = None
    
    ## Retrive the adata after data type selection
    if sample_name!= None and sample_name!= '---Please choose---':
        filename = sample_dic[sample_name]
        # adata = get_adata_aws(filename)
        # adata = get_adata_sftp(filename)
        adata = sc.read_h5ad('/Users/ziliangluo/Library/CloudStorage/OneDrive-UniversityofGeorgia/PycharmProjects/SpatialSeq/saved_ad/gma_sp_CS2A_fromSeurat.h5ad')
        
        # read the genes
        gene_ids = adata.var.index.tolist()
    else:
        st.sidebar.write('Please select a data to explore the genes')
    # print(filename)

    # gene_ids = ['test','ann1.Glyma.02G228100', 'ann1.Glyma.15G127900']
    
    if sample_name and sample_name != '---Please choose---':
        st.sidebar.markdown('## Please select gene to plot:')
        # plot_type = st.sidebar.selectbox('Select plot type', ['UMAP', 'Spatial'])
        plot_type = st.sidebar.radio('Select plot type', ['Violin','UMAP', 'Spatial'], horizontal=True)
        gene_name = st.sidebar.selectbox('Enter gene name for expression plot', ['', *gene_ids])
        
            
    ###############################################
    ##                Main page                 ###
    ###############################################
    st.markdown('## Soybean Gene Visualization')
    st.markdown('Welcome to the soybean multiomic single-cell database. Please select the dataset and genes from the sidebar to start. :balloon:')
    st.markdown('[A spatially resolved multiomic single-cell atlas of soybean development (paper link)](https://schmitzlab.uga.edu/), Zhang et al., 2024 BioRxiv')
    
    if sample_name and sample_name != '---Please choose---':
        # st.write('Generating the', plot_type, 'plot for', sample_name,  lib_type)
        
        # selected gene for plot
        variables_to_plot = ['cell_types']        ## cell typs not for Violin plot
        if gene_name:
            variables_to_plot.append(gene_name)
        
       
        # st.markdown('Selected gene `%s`' % gene_name)
        if plot_type == 'Violin':
            st.markdown('**Violin Plot**')
            st.markdown('Please select a gene.')
            if gene_name != '':
                fig1, ax1 = plt.subplots(figsize=(12, 6))
                sc.pl.violin(adata,gene_name, groupby='cell_types', rotation= 60, ax=ax1)
                st.pyplot(fig1)
            
        elif plot_type == 'UMAP':
            st.markdown('**UMAP Plot**')
            fig2, axs2 = plt.subplots(len(variables_to_plot), 1 , figsize=(5.5, 4.5* len(variables_to_plot)))
            if len(variables_to_plot) == 1:
                axs2 = [axs2]
            for ax, gene in zip(axs2, variables_to_plot):
                sc.pl.umap(adata, color=gene, ax=ax, show=False)
            plt.subplots_adjust(wspace=0.5)
            st.pyplot(fig2)
                        
        elif plot_type == 'Spatial':
            st.markdown('**Spatial Plot**')
            if lib_type != 'spRNA-seq':
                st.markdown(':exclamation: `Selected data is not spatial transcriptome. Please select other plot types`')
            else:
                fig3, axs3 = plt.subplots(len(variables_to_plot), 1 , figsize=(5.5, 4.5* len(variables_to_plot)))
                if len(variables_to_plot) == 1:
                    axs3 = [axs3]
                for ax, gene in zip(axs3, variables_to_plot):
                    sc.pl.spatial(adata, color=gene, ax=ax, show=False)
                plt.subplots_adjust(wspace=0.5)
                st.pyplot(fig3)
                    

if __name__ == '__main__':
    main()
