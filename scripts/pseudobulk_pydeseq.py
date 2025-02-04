import argparse
import scanpy as sc
import pandas as pd
from itertools import combinations
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
import decoupler as dc

def main():
    parser = argparse.ArgumentParser(description='Pseudobulk DESeq2 analysis for a specific tissue.')
    parser.add_argument('--tissue', '-t', type=str, help='Name of the tissue for analysis', required=True)
    parser.add_argument('--level', '-l', type=str, help='Annotation level: celltype, subtype, or general_celltype', required=True)
    
    args = parser.parse_args()
    
    tissue = args.tissue
    level = args.level
    
    file = f'../IGVF_analysis/cellbender_tissues/annotated/{tissue}_annotated.h5ad'
    
    adata = sc.read_h5ad(file)
    
    # Remove low quality cells
    adata = adata[adata.obs['celltype'] != 'low quality'].copy()
    
    # Count cells per cell type and genotype
    type_counts = adata.obs.groupby([level, 'Genotype']).size().unstack(fill_value=0)
    
    # Identify valid cell types
    valid_types = type_counts[type_counts.min(axis=1) >= 100].index.tolist()
    invalid_types = type_counts[type_counts.min(axis=1) < 100].index.tolist()
    
    print(f"Invalid cell types (less than 100 cells in any genotype): {invalid_types}")
    
    # Initialize list to collect results
    all_results = []
    
    # Loop through valid cell types
    for this_type in valid_types:
        print(f"Processing cell type: {this_type}")
        
        # Subset adata for the current cell type
        pdata_sub = adata[adata.obs[level] == this_type].copy()
        
        # Pseudobulking - decoupler with default options
        pbulk = dc.get_pseudobulk(pdata_sub,
                                  sample_col="lab_sample_id", 
                                  groups_col=None,
                                  layer='cellbender_counts',
                                  mode='sum',
                                  skip_checks=True,
                                  min_cells=10,
                                  min_counts=1000)
        
        # Save pseudobulk adata
        savetype = this_type.replace('+', '_').replace('/', '_').replace('-', '_').replace(' ', '_')
        output_dir = f'../IGVF_analysis/cellbender_tissues/pseudobulk/celltype_pseudobulk/{tissue}/'
        os.makedirs(output_dir, exist_ok=True)
        
        output_file = f'{tissue}_{savetype}_pseudobulk.h5ad'
        output_path = os.path.join(output_dir, output_file)
        pbulk.obs = pbulk.obs.applymap(str)
        pbulk.write_h5ad(output_path)
        
        # Get pbulk df
        df_pbulk = pbulk.to_df()
        
        # Extract metadata
        metadata = pbulk.obs[['lab_sample_id', 'Sex', 'Genotype']].copy()
        metadata = metadata.drop_duplicates()
        metadata.set_index('lab_sample_id', inplace=True)
        metadata = metadata.reindex(df_pbulk.index)
        
        # Get unique genotypes
        genotypes = pbulk.obs['Genotype'].unique().tolist()
        
        # Perform pairwise genotype comparisons
        for genotype1, genotype2 in combinations(genotypes, 2):
            print(f"Running DESeq2 for {this_type} comparing {genotype1} vs {genotype2}", flush=True)

            # Create the DESeq2 object and run the analysis
            if tissue in ["GonadsMale", "GonadsFemale"] or pbulk.obs["Sex"].nunique() < 2:
                dds = DeseqDataSet(counts=df_pbulk,
                                   metadata=metadata,
                                   design_factors=['Genotype'],
                                   ref_level=['Genotype', genotype2],
                                   refit_cooks=True, quiet=True)
            else:
                dds = DeseqDataSet(counts=df_pbulk,
                                   metadata=metadata,
                                   design_factors=['Sex', 'Genotype'],
                                   ref_level=['Genotype', genotype2],
                                   refit_cooks=True, quiet=True)

            # Compute LFCs
            dds.deseq2()

            stat_res = DeseqStats(dds, contrast=['Genotype', genotype1, genotype2])
            stat_res.summary()

            # Shrink LFCs
            coeff_name = f'Genotype_{genotype1}_vs_{genotype2}'
            stat_res.lfc_shrink(coeff=coeff_name)

            # Add comparison details to the results
            result_df = stat_res.results_df.reset_index()
            result_df['comparison'] = coeff_name
            result_df['celltype'] = this_type
            result_df['tissue'] = tissue
            
            print("Results:", result_df.head(), flush = True)

            # Append results to the list
            all_results.append(result_df)
    
    # Save all results to a CSV file
    output_dir = f'degs/{tissue}/'
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, 'results_de_analysis.csv')
    final_results = pd.concat(all_results, ignore_index=True)
    final_results.to_csv(output_file, index=False)
    print(f"All DESeq2 results saved to: {output_file}")

if __name__ == "__main__":
    main()

    
