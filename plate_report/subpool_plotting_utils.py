import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from natsort import natsorted


def plot_plate_map(df, color_by, color_dict, output_path="plots/plate_map_plot.png"):
    import matplotlib.pyplot as plt
    from matplotlib import patches
    import pandas as pd

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

        # Save mapping table
        multiplexed_key = pd.DataFrame(
            [(label, pair[0], pair[1]) for pair, label in pair_map.items()],
            columns=["Pair", "Multiplexed Genotype 1", "Multiplexed Genotype 2"]
        )
        multiplexed_key.to_csv("plots/multiplexed_genotype_key.csv", index=False)

    # Now plot using genotype
    rows = list('ABCDEFGH')
    columns = range(1, 13)
    color_map = {id: color_dict.get(id, 'white') for id in df[color_by].unique()}

    fig, ax = plt.subplots(figsize=(5, 3))

    for row_label in rows:
        for col_num in columns:
            well_data = df[(df['Row'] == row_label) & (df['Column'] == col_num)]
            if not well_data.empty:
                well_color = color_map[well_data[color_by].values[0]]
                sex = well_data['Sex'].values[0]
            else:
                well_color = 'white'
                sex = None
            ax.add_patch(plt.Circle((col_num, rows.index(row_label) + 1), 0.45, color=well_color, ec='black', linewidth=1))
            if sex == 'Female':
                ax.text(col_num, rows.index(row_label) + 1, 'X', color='black', ha='center', va='center', fontsize=8)

    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.5, 8.5)
    ax.set_xticks(columns)
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(rows, fontsize=11)
    ax.set_xticklabels(columns, fontsize=11)
    ax.set_title(f'Sample Loading Plate Map', fontsize=10)
    plt.gca().invert_yaxis()

    handles = [patches.Patch(color=color_dict[key], label=key) for key in color_dict]
    handles.append(patches.Patch(color='white', label='X = Female'))
    ax.legend(handles=handles, bbox_to_anchor=(1.05, 1.05), loc='upper left', fontsize=10, title=color_by)

    plt.tight_layout()
    fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    return output_path


def create_round_heatmap(obs, well_column, kit, sampletype):
    round_number = well_column[2]
    well_counts = obs[well_column].value_counts().sort_index()
    well_counts_df = well_counts.reset_index()
    well_counts_df.columns = ['well', 'cell_count']

    rows = list('ABCDEFGH')
    cols = list(range(1, 13))
    plate_df = pd.DataFrame(0, index=rows, columns=cols)

    for _, row in well_counts_df.iterrows():
        plate_df.at[row['well'][0], int(row['well'][1:])] = row['cell_count']

    if kit == "WT" and round_number == "1":
        plate_df = plate_df.loc['A':'D']

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.heatmap(
        plate_df, annot=True, fmt='d', cmap='Reds', cbar=True,
        linewidths=0, annot_kws={"size": 12},
        cbar_kws={'label': 'Cell Count'}, vmin=0, vmax=plate_df.max().max() * 1.2, ax=ax
    )
    ax.set_title(f'Round {round_number} - {sampletype} per well (>=500 UMI)', fontsize=16)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16, rotation=0)

    output_path = f"plots/round_{well_column}_plate_hmap_all_subs.png"
    fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    return output_path


def plot_knee_raw_counts(plate, subpool, output_path="plots/raw_counts_single_subpool_knee_plot_plate.png"):
    fig, ax = plt.subplots(figsize=(8, 6))

    file_path = f"/dfs9/seyedam-lab/erebboah/parse_pipeline/pipeline_report/{plate}/{subpool}/knee_plot_df.csv"
    knee_df = pd.read_csv(file_path)
    knee = knee_df['UMIs'].values
    ax.loglog(range(len(knee)), knee, linewidth=2, label=subpool)

    ax.set_ylabel("UMI Counts", fontsize=18)
    ax.set_xlabel("Barcodes", fontsize=18)
    ax.axhline(y=500, linewidth=2, color="#505050", linestyle='--', label='500 UMI')
    ax.axhline(y=200, linewidth=2, color="#A0A0A0", linestyle='--', label='200 UMI')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, which="both")

    fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    return output_path


def plot_knee_cb(obs, category_column='subpool', output_path="plots/cb_knee_subpool_plot.png"):
    subpool = obs[category_column].unique()[0]
    
    fig, ax = plt.subplots(figsize=(8, 6))

    knee_raw = np.sort(obs['total_counts_raw'].values)[::-1]
    knee_cb = np.sort(obs['total_counts_cb'].values)[::-1]
    ax.loglog(np.arange(len(knee_raw)), knee_raw, linewidth=3, label=f"{subpool} (raw)")
    ax.loglog(np.arange(len(knee_cb)), knee_cb, linewidth=3, label=f"{subpool} (cellbender)")

    ax.axhline(y=500, linewidth=2, color="#505050", linestyle='--', label='500 UMI')
    ax.axhline(y=200, linewidth=2, color="#A0A0A0", linestyle='--', label='200 UMI')
    ax.set_ylabel("UMI Counts", fontsize=18)
    ax.set_xlabel("Set of Barcodes", fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(fontsize=11, bbox_to_anchor=(1.02, 1.02), loc='upper left')
    ax.grid(True, which="both")

    fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    return output_path


def plot_stacked_barplots(obs, cluster_key, var_key, output_path, color_dict=None, figsize=(5, 12), annotations=True, reverse_order=False, cluster_order=None, title=None):
    grouped = obs.groupby([cluster_key, var_key]).size().unstack(fill_value=0)
    proportions = grouped.div(grouped.sum(axis=1), axis=0)

    if cluster_order:
        proportions = proportions.loc[cluster_order]
    elif reverse_order:
        proportions = proportions.iloc[::-1]

    cluster_sizes = grouped.sum(axis=1).loc[proportions.index]
    unique_cats = proportions.columns.tolist()
    colors = [color_dict.get(cat, "gray") for cat in unique_cats] if color_dict else sns.color_palette("husl", n_colors=len(unique_cats))

    fig, ax = plt.subplots(figsize=figsize)
    proportions.plot(kind="barh", stacked=True, color=colors, ax=ax, width=0.8)

    if annotations:
        for i, txt in enumerate(cluster_sizes):
            ax.text(1.02, i, f"{txt:,}", fontsize=12, va="center", transform=ax.get_yaxis_transform())

    ax.set_xlim(0, 1.1)
    ax.set_xlabel("Proportion", fontsize=14)
    ax.set_ylabel(cluster_key, fontsize=14)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_title(title if title else f"{var_key} by {cluster_key}", fontsize=16)
    if annotations:
        ax.legend(title=var_key, bbox_to_anchor=(1.2, 1), loc="upper left")
    else:
        ax.get_legend().remove()

    fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    return output_path


def plot_stacked_main(combined_obs, subpool, main_tissue, output_path):
    obs = combined_obs[combined_obs['Tissue'] == main_tissue]
    obs = obs[obs['subpool'] == subpool]
    return plot_stacked_barplots(
        obs=obs,
        cluster_key="Genotype",
        var_key="celltype",
        output_path=output_path,
        reverse_order=True,
        figsize=(6.5, 6),
        title=f"celltype by Genotype in {main_tissue}"
    )


def plot_stacked_mult(combined_obs, subpool, mult_tissue, output_path):
    obs = combined_obs[combined_obs['Tissue'] == mult_tissue]
    obs = obs[obs['subpool'] == subpool]
    return plot_stacked_barplots(
        obs=obs,
        cluster_key="Genotype",
        var_key="celltype",
        output_path=output_path,
        reverse_order=True,
        figsize=(6.5, 6),
        title=f"celltype by Genotype in {mult_tissue}"
    )


def plot_qc_violins(obs, output_prefix="plots/qc_violin_subpool_plot", size=14):
    if len(obs['Tissue'].unique()) > 5:
        plot_category = "plate"
    else:
        plot_category = "Tissue"

    # --- Plot 1: Number of genes and total counts ---
    data_ngb = pd.melt(obs, id_vars=[plot_category], 
                       value_vars=['n_genes_by_counts_raw', 'n_genes_by_counts_cb'], 
                       var_name='Count type', value_name='counts')
    data_ngb['Count type'] = data_ngb['Count type'].replace({
        'n_genes_by_counts_raw': 'Raw',
        'n_genes_by_counts_cb': 'CellBender'
    })

    data_tbc = pd.melt(obs, id_vars=[plot_category], 
                       value_vars=['total_counts_raw', 'total_counts_cb'], 
                       var_name='Count type', value_name='counts')
    data_tbc['Count type'] = data_tbc['Count type'].replace({
        'total_counts_raw': 'Raw',
        'total_counts_cb': 'CellBender'
    })

    fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharey=False)

    sns.violinplot(x=plot_category, y='counts', hue='Count type', data=data_ngb, split=True, ax=axes[0])
    axes[0].set_title('Number of Expressed Genes', fontsize=size)
    axes[0].set_ylabel('# genes', fontsize=size)
    axes[0].set_xlabel("", fontsize=size)
    axes[0].tick_params(axis='both', which='major', labelsize=size)

    sns.violinplot(x=plot_category, y='counts', hue='Count type', data=data_tbc, split=True, ax=axes[1])
    axes[1].set_title('Total Counts', fontsize=size)
    axes[1].set_xlabel("", fontsize=size)
    axes[1].set_ylabel('# Counts', fontsize=size)
    axes[1].tick_params(axis='both', which='major', labelsize=size)

    plt.tight_layout()
    fname1 = f"{output_prefix}1.png"
    plt.savefig(fname1, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)

    # --- Plot 2: Mitochondrial percent and doublet score ---
    data_pmt = pd.melt(obs, id_vars=[plot_category], 
                       value_vars=['pct_counts_mt_raw', 'pct_counts_mt_cb'], 
                       var_name='Count type', value_name='Percent')
    data_pmt['Count type'] = data_pmt['Count type'].replace({
        'pct_counts_mt_raw': 'Raw',
        'pct_counts_mt_cb': 'CellBender'
    })

    data_ds = pd.melt(obs, id_vars=[plot_category], 
                      value_vars=['doublet_score'], 
                      var_name='Count type', value_name='Score')
    data_ds['Count type'] = 'Raw'

    fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharey=False)

    sns.violinplot(x=plot_category, y='Percent', hue='Count type', data=data_pmt, split=True, ax=axes[0])
    axes[0].set_title('% Mitochondrial Expression', fontsize=size)
    axes[0].set_ylabel('Percent', fontsize=size)
    axes[0].set_xlabel("", fontsize=size)
    axes[0].tick_params(axis='both', which='major', labelsize=size)

    sns.violinplot(x=plot_category, y='Score', hue='Count type', data=data_ds, ax=axes[1])
    axes[1].set_title('Scrublet Doublet Score', fontsize=size)
    axes[1].set_ylabel('Score', fontsize=size)
    axes[1].set_xlabel("", fontsize=size)
    axes[1].tick_params(axis='both', which='major', labelsize=size)

    plt.tight_layout()
    fname2 = f"{output_prefix}2.png"
    plt.savefig(fname2, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)

    return fname1, fname2


def plot_qc_violins_filtered(obs_filt, output_prefix="plots/qc_violin_subpool_plot_filt", size=14):
    return plot_qc_violins(obs_filt, output_prefix=output_prefix, size=size)
