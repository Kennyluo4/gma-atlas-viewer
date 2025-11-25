import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import boto3        # for s3 file transfer
# import paramiko     # for sftp file transfer
# from io import StringIO
import os

# Set the page config to wide mode
st.set_page_config(layout="wide")
    
# Inject custom CSS to set specific width or ratio
st.markdown(
    """
    <style>
    /* Define the main content width */
    .main {
        max-width: 80%; /* Adjust the percentage as needed */
        margin: 1 auto; /* Center align */
    }
    /* Optional: Define sidebar width */
    .sidebar .sidebar-content {
        width: 20%; /* Adjust the percentage as needed */
    }
    </style>
    """,
    unsafe_allow_html=True
)

@st.cache_data
def get_adata_aws(file_name):
    '''read files stored on AWS S3'''
    # Read AWS credentials from environment variables
    # aws_access_key_id = st.secrets['AWS_ACCESS_KEY_ID']           # st.secrets works only for streamlit cloud hosted app
    # aws_secret_access_key = st.secrets['AWS_SECRET_ACCESS_KEY']
    aws_access_key_id = os.getenv('AWS_ACCESS_KEY_ID')
    aws_secret_access_key = os.getenv('AWS_SECRET_ACCESS_KEY')

    # Initialize S3 client
    s3 = boto3.client('s3', aws_access_key_id=aws_access_key_id, aws_secret_access_key=aws_secret_access_key)

    bucket_name = 'soybean-atlas'
    file_key = file_name
    local_file_name = file_name
    print(f'Reading file: {local_file_name}')

    # Download file from Se`3
    s3.download_file(bucket_name, file_key, local_file_name)

    # Read the AnnData file using Scanpy
    adata = sc.read(local_file_name)
    return adata

# def get_adata_sftp(file_name):
#     '''read files stored on cluster'''
#     hostname = 'xfer.gacrc.uga.edu'
#     # st.secrets['xferhost']
#     port = 22
#     username = ''
#     # st.secrets['userid']
#     password = ''
#     # st.secrets['xferpass']
#     remote_path = '/' + file_name
#     with open(os.path.expanduser('~/.ssh/id_rsa'), 'r') as file:
#         key_content = file.read()
#     try:
#         private_key = paramiko.RSAKey.from_private_key(StringIO(key_content), password='')

#         transport = paramiko.Transport((hostname, port))
#         transport.connect(username=username, pkey=private_key)

#         sftp = paramiko.SFTPClient.from_transport(transport)
#         with sftp.file(remote_path, 'r') as remote_file:
#             file_content = remote_file.read()
#         sftp.close()
#         transport.close()
#         return file_content
#     except Exception as e:
#         st.error(f"An error occurred: {e}")
#         return None

def main():
    ## dic for matching input sample same and stored adata file name
    sample_dic = {'heart stage seed (spatial)': 'sp_adt_HS_A_slim_062424.h5ad','cotyledon stage seed (spatial)':'sp_adt_CS_A_slim_062424.h5ad', 
                  'early maturation stage seed (spatial)':'sp_adt_ES_A_slim_062424.h5ad','hypocotyl (spatial)': 'sp_adt_HP_A_slim_062424.h5ad','root (spatial)': 'sp_adt_RTB_rep1_slim_062424.h5ad',
                  
                  'globular stage seed (snRNA)': 'snrna_adt_GS_slim.h5ad','heart stage seed (snRNA)': 'snrna_adt_HS_slim.h5ad',
                  'cotyledon stage seed (snRNA)': 'snrna_adt_CS_slim.h5ad', 'early maturation stage seed (snRNA)': 'snrna_adt_ES_slim.h5ad',
                  'hypocotyl (snRNA)':'snrna_adt_HP_slim.h5ad', 'root (snRNA)':'snrna_adt_RT_slim.h5ad', 'nodule (snRNA)':'snrna_adt_NOD_slim.h5ad',
                  'root with nematode infection (snRNA)':'Gm_All_root_integrated-infected_RNA_RNA_imputed.h5ad','root with mock infection (snRNA)':'Gm_All_root_integrated-mock_RNA_RNA_imputed.h5ad','All roots nematode infected and mock (snRNA)':'Gm_All_root_integrated_RNA_RNA_imputed.h5ad',
                  
                  'globular stage seed (motif)' : 'gma_GS_motif_deviation_imputed_slim.h5ad','heart stage seed (motif)': 'gma_HS_motif_deviation_imputed_slim.h5ad',
                  'cotyledon stage seed (motif)': 'gma_CS_motif_deviation_imputed_slim.h5ad', 'early maturation stage seed (motif)': 'gma_ES_motif_deviation_imputed_slim.h5ad', 'middle maturation stage seed (motif)':'gma_MS_motif_deviation_imputed_slim.h5ad',
                  'hypocotyl (motif)':'gma_HP_motif_deviation_imputed_slim.h5ad', 'root (motif)':'gma_RT_motif_deviation_imputed_slim.h5ad', 'nodule (motif)':'gma_NOD_motif_deviation_imputed_slim.h5ad',
                  'leaf (motif)':'gma_LF_motif_deviation_imputed_slim.h5ad', 'pod (motif)':'gma_PD_motif_deviation_imputed_slim.h5ad',
                  
                  'globular stage seed (TF)': 'gma_GS_TF_expression_imputed_slim.h5ad','heart stage seed (TF)': 'gma_HS_TF_expression_imputed_slim.h5ad',
                  'cotyledon stage seed (TF)': 'gma_CS_TF_expression_imputed_slim.h5ad', 'early maturation stage seed (TF)': 'gma_ES_TF_expression_imputed_slim.h5ad',
                  'hypocotyl (TF)':'gma_HP_TF_expression_imputed_slim.h5ad', 'root (TF)':'gma_RT_TF_expression_imputed_slim.h5ad', 'nodule (TF)':'gma_NOD_TF_expression_imputed_slim.h5ad',
                  }
    
    ###############################################
    ##                 Sidebar                  ###
    ###############################################
    # with st.sidebar.form():           # use if want to add a Submit button before change everything
    
    st.sidebar.header('Plot Configuration')
    st.sidebar.markdown('## Please select a dataset:')
    lib_type = st.sidebar.selectbox('Data', [ '---Please choose---','snRNA-seq', 'spRNA-seq', 'scATAC-seq'])
    
    if lib_type == 'snRNA-seq':
        sample_name = st.sidebar.selectbox('Tissue', ['---Please choose---', 'globular stage seed','heart stage seed', 'cotyledon stage seed',
                                                      'early maturation stage seed','hypocotyl','root','nodule', 'root with nematode infection','root with mock infection','All roots nematode infected and mock'])
        sample_name = sample_name + ' (snRNA)'
    elif lib_type == 'spRNA-seq':
        sample_name = st.sidebar.selectbox('Tissue', ['---Please choose---', 'heart stage seed', 'cotyledon stage seed', 'early maturation stage seed',
                                                      'hypocotyl','root'])
        sample_name = sample_name + ' (spatial)'
    elif lib_type == 'scATAC-seq':
        data_type = st.sidebar.radio('Select data type', ['motif deviation', 'motif associated TF expression'], horizontal=True)
        if data_type =='motif deviation':
            sample_name = st.sidebar.selectbox('Tissue', ['---Please choose---', 'globular stage seed','heart stage seed', 'cotyledon stage seed',
                                                        'early maturation stage seed','middle maturation stage seed','hypocotyl','root','nodule',
                                                        'leaf', 'pod'])
            sample_name = sample_name + ' (motif)'
        elif data_type =='motif associated TF expression':
            sample_name = st.sidebar.selectbox('Tissue', ['---Please choose---', 'globular stage seed (TF)','heart stage seed', 'cotyledon stage seed',
                                                        'early maturation stage seed','hypocotyl','root','nodule'])
            sample_name = sample_name + ' (TF)'
    else:
        sample_name = None
    
    ## Retrive the adata after data type selection
    if sample_name!= None and sample_name!= '---Please choose---':
        filename = sample_dic[sample_name]
        adata = sc.read_h5ad(filename)
        # adata = get_adata_aws(filename)
        # adata = get_adata_sftp(filename)
        # adata = sc.read_h5ad('/Users/ziliangluo/Library/CloudStorage/OneDrive-UniversityofGeorgia/PycharmProjects/SpatialSeq/saved_ad/Gm_atlas_Cotyledon_stage_seeds_ATAC_imputed2.h5ad')
        
        # read the genes
        gene_ids = adata.var.index.tolist()
    else:
        st.sidebar.write('Please select a data to explore the genes')
    # print(filename)

    # gene_ids = ['test','ann1.Glyma.02G228100', 'ann1.Glyma.15G127900']
    
    plot_type = None
    if sample_name and sample_name != '---Please choose---':
        st.sidebar.markdown('## Please select gene to plot:')
        # plot_type = st.sidebar.selectbox('Select plot type', ['UMAP', 'Spatial'])
        plot_type = st.sidebar.radio('Select plot type', ['UMAP', 'Spatial','Violin'], horizontal=True)
        gene_name = st.sidebar.selectbox('Enter gene name for expression plot', ['', *gene_ids])
        
            
    ###############################################
    ##                Main page                 ###
    ###############################################
    
    st.markdown('## Soybean Gene Viewer')
    
    if sample_name and sample_name != '---Please choose---':
        # st.write('Generating the', plot_type, 'plot for', sample_name,  lib_type)
        
        # selected gene for plot
        variables_to_plot = ['cell_types']        ## cell typs not for Violin plot
        if gene_name:
            variables_to_plot.append(gene_name)
        # st.markdown('Selected gene `%s`' % gene_name)
        
        if plot_type == 'UMAP':
            st.subheader('UMAP Plot')
            fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(5, 5 * len(variables_to_plot)))
            if len(variables_to_plot) == 1:
                axs = [axs]
            for ax, gene in zip(axs, variables_to_plot):
                sc.pl.umap(adata, color=gene, ax=ax, show=False)
                plt.subplots_adjust(wspace=1.2)
            st.pyplot(fig)
        elif plot_type == 'Spatial':
            st.subheader('Spatial Plot')
            if lib_type != 'spRNA-seq':
                st.error(' :crying_cat_face: Spatial Map not available. \nSelected data is not spatial transcriptomics')
            else:
                fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(5, 5 * len(variables_to_plot)))
                if len(variables_to_plot) == 1:
                    axs = [axs]
                for ax, gene in zip(axs, variables_to_plot):
                    sc.pl.spatial(adata, color=gene, ax=ax, show=False)
                    plt.subplots_adjust(wspace=1.2)
                st.pyplot(fig)
        elif plot_type == 'Violin':
            st.markdown('**Violin Plot**')
            st.error(':point_left: Please select a gene in the sidebar')
            if gene_name != '':
                fig, ax = plt.subplots(figsize=(12, 6))
                sc.pl.violin(adata,gene_name, groupby='cell_types', rotation= 60, ax=ax)
                # Add custom y-axis label
                ax.set_ylabel(f'{gene_name} (normalized expression)')  
                st.pyplot(fig)
                
    if plot_type == None:
        st.markdown('Welcome to the soybean multi-omic single-cell database. Please select the dataset and genes from the sidebar to start.')
        st.markdown('[A spatially resolved multi-omic single-cell atlas of soybean development](https://doi.org/10.1016/j.cell.2024.10.050), Zhang et al., 2024 Cell')
        st.image('./images/overview.png')
        st.image('./images/nematode.png')
    ## footnote
    footnote = """
                <hr>
                <p style='font-size: small;'>For issues and questions, please contact 
                <a href='mailto:luoziliang@uga.edu'>Ziliang</a>.</p>
                """
    st.markdown(footnote, unsafe_allow_html=True)
    
    # if sample_name == None:
    #     sentiment_mapping = ["one", "two", "three", "four", "five"]
    #     st.markdown('Rate our website')
    #     selected = st.feedback('stars')
    #     if selected is not None:
    #         st.markdown(f'You selected {sentiment_mapping[selected]} star(s).')
    #         record = str(selected) + '\n'
    #         with open('feedbacks.txt', 'a') as handle:
    #             handle.write(record)
                
if __name__ == '__main__':
    main()
