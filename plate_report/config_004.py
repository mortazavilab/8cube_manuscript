# ----- Manual Inputs -----
plate = 'igvf_004'
kit = 'WT_mega'
chemistry = 'v2'
sampletype = 'nuclei'

min_counts = 500
max_counts = 150000
min_genes = 250
pct_counts_mt = 1
doublet_score = 0.25

# ----- File Paths -----
config_path = f'/dfs9/seyedam-lab/erebboah/parse_pipeline/configs/{plate}_config_with_exome.tsv'
sample_metadata_path = '/dfs9/seyedam-lab/erebboah/parse_pipeline/configs/sample_metadata.csv'
obs_path = f'{plate}/{plate}_adata_obs.csv'