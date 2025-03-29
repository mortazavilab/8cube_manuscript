import pandas as pd
import numpy as np
import re
import os


def parse_sample_df(fname):
    df = pd.read_csv(fname)
    df.rename({'Experiment': 'plate'}, axis=1, inplace=True)
    g_cols = ['mult_genotype_1', 'mult_genotype_2']
    df[g_cols] = df.Genotype.str.split('/', expand=True)
    df.loc[df.well_type == 'Single', g_cols] = np.nan
    inds = df.loc[df.well_type == 'Multiplexed'].index
    df['mult_genotype'] = np.nan
    df.loc[inds, 'mult_genotype'] = df.loc[inds, g_cols[0]].astype(str) + '_' + df.loc[inds, g_cols[1]].astype(str)
    return df


def get_subpools(config_df, plate):
    subpools = config_df['subpool'].unique().tolist()
    regular_subpools = sorted(
        [s for s in subpools if s != 'Subpool_EX'],
        key=lambda x: int(re.search(r'\d+', x).group())
    )

    # Define sets of plates for conditions
    add_ex_plates = {'igvf_003', 'igvf_004', 'igvf_005', 'igvf_007', 'igvf_008b', 'igvf_009', 'igvf_010', 'igvf_011', 'igvf_012'}
    regular_only_plates = {'igvf_008', 'igvf_012', 'igvf_015'}

    if plate in regular_only_plates:
        return regular_subpools
    elif plate in add_ex_plates:
        return ['Subpool_EX'] + regular_subpools
    else:
        return regular_subpools


def get_color_dict(df, color_by):
    # Default to Genotype
    df['genotype'] = df['Genotype']

    # If multiplexed wells exist...
    if 'well_type' in df.columns and any(df['well_type'] == 'Multiplexed'):
        pair_cols = ['mult_genotype_1', 'mult_genotype_2']
        multiplexed = df[df['well_type'] == 'Multiplexed']

        # Get sorted unique genotype pairs
        multiplexed_pairs = (
            multiplexed[pair_cols]
            .apply(lambda row: tuple(sorted(row)), axis=1)
            .drop_duplicates()
            .reset_index(drop=True)
        )

        # Map to letters A, B, C...
        pair_labels = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")[:len(multiplexed_pairs)]
        pair_map = {pair: label for pair, label in zip(multiplexed_pairs, pair_labels)}

        # Assign label to each row
        def assign_label(row):
            if row['well_type'] == 'Multiplexed':
                pair = tuple(sorted((row['mult_genotype_1'], row['mult_genotype_2'])))
                return pair_map.get(pair, 'Unknown')
            return row['Genotype']

        df['genotype'] = df.apply(assign_label, axis=1)
        
    geno_dict = {
        'B6J': '#C0BFBF', 
        'NODJ': '#4F6EAF', 
        'AJ': '#F4C245', 
        'PWKJ': '#D83026', 
        '129S1J': '#DA9CC1',
        'CASTJ': '#55AF5B', 
        'WSBJ': '#683C91', 
        'NZOJ': '#52A5DB', 
        'A': '#1b9f78',
        'B': '#da5f00', 
        'C': '#766fb3', 
        'D': '#e8298a', 
        'E': '#67a51e',
        'F': '#e6ab01', 
        'G': '#9B958D', 
        'H': '#1f77b5',
        'B6CASTF1J': '#7e9772',
        'B6PWKF1J': '#a5655f', 
        'B6WSBF1J': '#827298', 
        'B6NODF1J': '#7880a1', 
        'B6AF1J': '#bea468',
        'B6NZOF1J': '#a7761d', 
        'B6129S1F1J': '#b893a8', 

    }
    tissue_dict = {
        'CortexHippocampus': '#F8C471', 
        'Heart': '#AF7AC5', 
        'Liver': '#F15922', 
        'DiencephalonPituitary': '#27AE60',
        'Adrenal': '#8B4513', 
        'Kidney': '#CF1E1A', 
        'Gastrocnemius': '#FF69B4', 
        'GonadsFemale': '#2FE3E5', 
        'GonadsMale': '#00B7B8'
    }
    dict_to_use = geno_dict if color_by == 'genotype' else tissue_dict
    values = df[color_by].dropna().unique()
    return {val: dict_to_use.get(val, 'white') for val in values}


def summarize_cell_counts(qc_200, qc_500, align_df, sampletype):
    n_cells_col = f'Number of {sampletype}'
    cells_200 = sum(int(x.replace(',', '')) for x in qc_200[n_cells_col])
    cells_500 = sum(int(x.replace(',', '')) for x in qc_500[n_cells_col])
    reads = sum(int(x.replace(',', '')) for x in align_df['Total reads'])
    return (f"{cells_200:,}", f"{cells_500:,}", f"{reads:,}")


def load_alignment_stats(plate, subpool):
    file_path = f"/dfs9/seyedam-lab/erebboah/parse_pipeline/pipeline_report/{plate}/{subpool}/report.csv"
    df = pd.read_csv(file_path)
    df['Subpool'] = subpool
    df_selected = df[[
        'Subpool', 'Total reads', 'Total UMIs', 'Sequencing saturation',
        'Total pseudoaligned', 'Pct. pseudoaligned']].copy()
    df_selected = df_selected.rename(columns={'Sequencing saturation': 'Fraction duplicated UMIs'})
    return df_selected


def load_qc_stats(plate, subpool, sampletype):
    def load_and_rename(file_path, umi_thresh):
        df = pd.read_csv(file_path)
        df['Subpool'] = subpool
        subset = df[[
            f'Subpool',
            f'Num. cells >= {umi_thresh} UMI',
            f'Median UMIs/cell ({umi_thresh} UMIs)',
            f'Mean UMIs/cell ({umi_thresh} UMIs)',
            f'Median genes/cell ({umi_thresh} UMIs)',
            f'Mean genes/cell ({umi_thresh} UMIs)']].copy()
        return subset.rename(columns={
            f'Num. cells >= {umi_thresh} UMI': f'Number of {sampletype}',
            f'Median UMIs/cell ({umi_thresh} UMIs)': 'Median UMIs',
            f'Mean UMIs/cell ({umi_thresh} UMIs)': 'Mean UMIs',
            f'Median genes/cell ({umi_thresh} UMIs)': 'Median genes',
            f'Mean genes/cell ({umi_thresh} UMIs)': 'Mean genes'
        })
    file_path = f"/dfs9/seyedam-lab/erebboah/parse_pipeline/pipeline_report/{plate}/{subpool}/report.csv"    
    qc_200 = load_and_rename(file_path, 200)
    qc_500 = load_and_rename(file_path, 500)
    return qc_200, qc_500


def filter_obs(obs, min_counts, max_counts, min_genes, pct_counts_mt, doublet_score):
    return obs[
        (obs['total_counts_raw'] > min_counts) &
        (obs['total_counts_raw'] < max_counts) &
        (obs['n_genes_by_counts_cb'] > min_genes) &
        (obs['pct_counts_mt_cb'] < pct_counts_mt) &
        (obs['doublet_score'] < doublet_score)
    ]


def extract_main_tissue(obs, plate):
    return obs[(obs['plate'] == plate) & (obs['well_type'] == 'Single')]['Tissue'].iloc[0]


def extract_mult_tissue(obs, plate):
    return obs[(obs['plate'] == plate) & (obs['well_type'] == 'Multiplexed')]['Tissue'].iloc[0]


def load_combined_obs(sample_df):
    tissues = sample_df['Tissue'].unique().tolist()
    all_obs = []
    
    for tissue in tissues:
        filepath = f"/share/crsp/lab/seyedam/erebboah/parse_pipeline/IGVF_analysis/cellbender_tissues_exome/obs_tables/{tissue}_integrated_processed_annotated_metadata.csv"
        df = pd.read_csv(filepath)
        df['Tissue'] = tissue
        all_obs.append(df)
    return pd.concat(all_obs, axis=0, ignore_index=True)


def load_cellbender_stats(plate, subpool, sampletype):
    all_cb_results = []
    metrics_path = f"/dfs10/bio/erebboah/parse_pipeline/pipeline_output_cellbender/{plate}/{subpool}/adata_denoised_metrics.csv"
    cb_metrics = pd.read_csv(metrics_path, header=None)

    selected_metrics = cb_metrics[cb_metrics[0].isin([
        'fraction_counts_removed_from_cells', 'found_cells', 'output_average_counts_per_cell'
    ])].copy()

    selected_metrics.loc[selected_metrics[0] == 'fraction_counts_removed_from_cells', 1] *= 100
    selected_metrics[0] = selected_metrics[0].replace({
        'fraction_counts_removed_from_cells': 'Percent counts removed',
        'found_cells': 'Found cells',
        'output_average_counts_per_cell': 'Mean denoised counts'
    })

    selected_metrics.loc[selected_metrics[0] == 'Found cells', 1] = selected_metrics.loc[selected_metrics[0] == 'Found cells', 1].apply(lambda x: f"{int(x):,}")
    selected_metrics.loc[selected_metrics[0] == 'Mean denoised counts', 1] = selected_metrics.loc[selected_metrics[0] == 'Mean denoised counts', 1].apply(lambda x: f"{x:,.1f}")
    selected_metrics.loc[selected_metrics[0] == 'Percent counts removed', 1] = selected_metrics.loc[selected_metrics[0] == 'Percent counts removed', 1].apply(lambda x: f"{x:,.1f}")

    cb_result = selected_metrics.set_index(0).T
    cb_result['Subpool'] = subpool

    cb_result = cb_result[['Subpool', 'Percent counts removed', 'Found cells', 'Mean denoised counts']]
    cb_result = cb_result.rename(columns={'Found cells': f'Found {sampletype}'})
    return cb_result

def load_cellbender_settings(plate, subpool):

    path = f"/dfs9/seyedam-lab/erebboah/parse_pipeline/pipeline_report/{plate}/{subpool}/cellbender_settings.csv"
    cb_settings = pd.read_csv(path)
    cb_settings = cb_settings.drop(columns=['plate'], errors='ignore')
    cb_settings = cb_settings.fillna("NA")
    cb_settings = cb_settings[cb_settings['subpool'] == subpool]
    cb_settings = cb_settings.rename(columns={'subpool': 'Subpool'})
    return cb_settings

def create_barcode_sample_map(plate, kit, chemistry):
    plate_upper = plate.upper()
    acc_map = pd.read_csv(f"/dfs9/seyedam-lab/erebboah/parse_pipeline/configs/{plate_upper}_barcode_sample_map.tsv", sep='\t')

    bcs = pd.read_csv(f'/dfs9/seyedam-lab/erebboah/parse_pipeline/ref/r1_RT_replace_{kit}_{chemistry}.txt', sep=' ', header=None)
    bcs.columns = ['R', 'T']
    bcs['T'] = bcs['T'].str.replace('*', '', regex=False)
    bcs = bcs.melt(value_vars=['R', 'T'], var_name='parse barcode type', value_name='barcode')
    acc_map = acc_map.merge(bcs, on='barcode', how='left')

    acc_map = acc_map.rename(columns={'accession': 'sample accession', 
                                      'tissue': 'sample description'})
    acc_map['position'] = 17
    acc_map = acc_map[['barcode', 'sample accession', 'sample description', 'position', 'parse barcode type', 'well']]
    acc_map.to_csv(f'{plate}/{plate}_barcode_sample_map.csv', index = False)
    #acc_map['sample accession'] = acc_map['sample accession'].str.replace(',', '\n')
    return acc_map