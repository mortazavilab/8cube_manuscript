import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import random
import anndata
import scipy as sp
import argparse
import os

# first run
#python3 format_adata_for_topyfic.py -d ../IGVF_analysis/cellbender_tissues/annotated/Satellite.h5ad -t Satellite

def main():
    parser = argparse.ArgumentParser(description='Regulatory topics modeling using Topyfic')
    parser.add_argument('--data', '-d', type=str, help='Path to input adata.', required=True)
    parser.add_argument('--tissue', '-t', type=str, help='Tissue name, first letter capital, no abbreviations; e.g. Gastrocnemius', required=True)
    parser.add_argument('--k', '-k', type=int, help='k value for Topyfic; higher number means more topics.', required=True)
    parser.add_argument('--random_state', '-r', type=int, help='Random state initialization.', required=True)
    
    args = parser.parse_args()
    
    path = args.data
    tissue = args.tissue
    k = args.k
    random_state = args.random_state
    
    # read in annotated adata
    data = sc.read_h5ad(path)

    # run topyfic
    train = Topyfic.Train(name=f'{tissue}_{k}_{random_state}',
                          k=k,
                          n_runs=1,
                          random_state_range=[random_state])
    

    print("Training model...")
    train.run_LDA_models(data, n_jobs=1, n_thread=1)
    
    if not os.path.exists(tissue):
        os.makedirs(tissue)
    
    train.save_train(save_path = f'{tissue}/')
    
if __name__ == "__main__":
    main()    
