import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import random
import anndata
import scipy as sp
import argparse

def normalization(mtx):
    pf = mtx.sum(axis=1).A.ravel()
    log1p_pf = np.log1p(sp.sparse.diags(pf.mean()/pf) @ mtx)
    
    pf = log1p_pf.sum(axis=1).A.ravel()
    pf_log1p_pf = sp.sparse.diags(pf.mean()/pf) @ log1p_pf
    
    return pf_log1p_pf


def main():
    
    parser = argparse.ArgumentParser(description='Prepare input adata for Topyfic; subset for regulatory genes, normalize, and save.')
    parser.add_argument('--data', '-d', type=str, help='Path to annotated adata.', required=True)
    parser.add_argument('--tissue', '-t', type=str, help='Tissue name, first letter capital, no abbreviations; e.g. Gastrocnemius', required=True)
  
    args = parser.parse_args()
    
    path = args.data
    tissue = args.tissue

    # read in annotated adata
    adata = sc.read_h5ad(path)

    # subset for regulatory genes
    adata_var = adata.var
    reg_genes = pd.read_csv("/share/crsp/lab/seyedam/share/igvf_pipeline/topyfic/gene_display_table_7categories.csv", index_col = 0)
    reg_genes['regulatory_genes'] = True
    reg_genes = reg_genes[["regulatory_genes", "biotype", "gene_name"]]

    adata_var = adata_var.reset_index(drop=True)
    merged_df = pd.merge(adata_var, reg_genes, how='left', left_on='gene_name', right_on='gene_name')
    merged_df['regulatory_genes'].fillna(False, inplace=True)
    merged_df.index = adata.var.index

    adata.var['regulatory_genes'] = merged_df['regulatory_genes']   
    adatas = adata[:, adata.var["regulatory_genes"]]

    # normalize
    X = adatas.layers['cellbender_counts'] # CELLBENDER COUNTS!
    data = anndata.AnnData(X=X, obs=adatas.obs, var=adatas.var)    
    data.var.set_index('gene_name', inplace=True)
    data.X = normalization(data.X)
    data.X = round(data.X)

    print(data)
    
    data.write_h5ad(f'{tissue}_regulatory_only_normalized.h5ad')
    
    
if __name__ == "__main__":
    main()    
