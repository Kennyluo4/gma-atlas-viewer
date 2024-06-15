import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import boto3

@st.cache_data
def get_adata(file_name):
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
    
    ## Retrive the adata after selection
    ##get the adata from s3
    try:
        filename = sample_dic[sample_name]
    except KeyError:
        st.write('Please select the data')
    # print(filename)
    adata = get_adata(filename)
    # adata = sc.read_h5ad('/Users/ziliangluo/Library/CloudStorage/OneDrive-UniversityofGeorgia/PycharmProjects/SpatialSeq/saved_ad/gma_sp_CS2A_fromSeurat.h5ad')
     
    gene_ids = adata.var.index.tolist()
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
            if gene_name != '':
                fig, ax = plt.subplots(figsize=(12, 6))
                sc.pl.violin(adata,gene_name, groupby='cell_types', rotation= 60, ax=ax)
                st.pyplot(fig)
            
        elif plot_type == 'UMAP':
            st.markdown('**UMAP Plot**')
            fig, axs = plt.subplots(len(variables_to_plot), 1 , figsize=(5.5, 4.5* len(variables_to_plot)))
            if len(variables_to_plot) == 1:
                axs = [axs]
            for ax, gene in zip(axs, variables_to_plot):
                sc.pl.umap(adata, color=gene, ax=ax, show=False)
            # plt.subplots_adjust(wspace=0.5)
            st.pyplot(fig)
                        
        elif plot_type == 'Spatial':
            st.markdown('**Spatial Plot**')
            if lib_type != 'spRNA-seq':
                st.markdown(':exclamation: `Selected data is not spatial transcriptome. Please select other plot types`')
            else:
                fig, axs = plt.subplots(len(variables_to_plot), 1 , figsize=(5.5, 4.5* len(variables_to_plot)))
                if len(variables_to_plot) == 1:
                    axs = [axs]
                for ax, gene in zip(axs, variables_to_plot):
                    sc.pl.spatial(adata, color=gene, ax=ax, show=False)
                # plt.subplots_adjust(wspace=0.5)
                st.pyplot(fig)
                    

if __name__ == '__main__':
    main()
