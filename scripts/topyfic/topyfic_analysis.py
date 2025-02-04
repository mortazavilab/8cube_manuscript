import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import random
import anndata
import scipy as sp
import argparse

def main():
    parser = argparse.ArgumentParser(description='Regulatory topics modeling using Topyfic')
    parser.add_argument('--data', '-d', type=str, help='Path to input adata.', required=True)
    parser.add_argument('--tissue', '-t', type=str, help='Tissue name, first letter capital, no abbreviations; e.g. Gastrocnemius', required=True)
    parser.add_argument('--k', '-k', type=int, help='k value for Topyfic; higher number means more topics.', required=True)
    
    args = parser.parse_args()
    
    path = args.data
    tissue = args.tissue
    k = args.k

    print("Combining train objects...")
    data = sc.read_h5ad(path)
    X = data.X
    obs = data.obs
    var = data.var

    adata = anndata.AnnData(X=X, obs=obs, var=var) 
    main_train = Topyfic.Train(name=f"{tissue}_" + str(k),
                           k=k,
                           n_runs=100)
    trains = []
    for i in range(100):
        train = Topyfic.read_train(f"{tissue}/train_{tissue}_{k}_{i}.p")
        trains.append(train)
    main_train.combine_LDA_models(data, single_trains=trains)
    main_train.save_train(save_path=f'{tissue}/')

    print("Creating TopModel...")
    train = Topyfic.read_train(f"{tissue}/train_{tissue}_{k}.p")
    top_model, clustering, adata = Topyfic.calculate_leiden_clustering(trains=[train], data=adata)
    top_model.save_topModel(save_path=f"{tissue}/")

    print("Analyzing TopModel...")    
    top_model = Topyfic.read_topModel(f"{tissue}/topModel_{tissue}_{k}.p")
    analysis_top_model = Topyfic.Analysis(Top_model=top_model)
    analysis_top_model.calculate_cell_participation(data=data)
    analysis_top_model.save_analysis(f"{tissue}/analysis_{tissue}_{k}")


    
if __name__ == "__main__":
    main()    