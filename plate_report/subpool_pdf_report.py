from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
import pandas as pd

def create_pdf_doc(plate, subpool):
    pdf_filename = f"{plate}/{plate}_{subpool}_report.pdf"
    return SimpleDocTemplate(
        pdf_filename,
        pagesize=letter,
        title=f"{plate.upper()} {subpool} Pipeline Report"
    )

def add_title(elements, title):
    title_style = ParagraphStyle(name='Title', fontSize=20, spaceAfter=20)
    elements.append(Paragraph(title, title_style))
    elements.append(Spacer(1, 5))

def add_section_header(elements, title):
    subtitle_style = ParagraphStyle(name='Subtitle', fontSize=16, spaceAfter=10)
    elements.append(Paragraph(title, subtitle_style))
    elements.append(Spacer(1, 10))

def add_subtitle(elements, subtitle):
    subsubtitle_style = ParagraphStyle(name='Subsubtitle', fontSize=12, textColor=HexColor('#666666'), spaceAfter=10)
    elements.append(Paragraph(subtitle, subsubtitle_style))
    elements.append(Spacer(1, 5))

def add_table(elements, df, title=None):
    if title:
        add_subtitle(elements, title)
    data = [df.columns.tolist()] + df.values.tolist()
    table = Table(data)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))
    elements.append(table)
    elements.append(Spacer(1, 10))

def add_image(elements, img_path, width, height):
    img = Image(img_path)
    img.drawWidth = width
    img.drawHeight = height
    img.hAlign = 'CENTER'
    elements.append(img)
    elements.append(Spacer(1, 10))

def add_paragraph(elements, text, style_name='BodyText', font_size=12, space_after=10, leading=15):
    style = ParagraphStyle(name=style_name, fontSize=font_size, spaceAfter=space_after, leading=leading)
    elements.append(Paragraph(text, style))
    elements.append(Spacer(1, 5))

def add_layout(elements, left_table, right_image_path):
    image = Image(right_image_path)
    image.drawWidth = 292.5
    image.drawHeight = 225
    image.hAlign = 'LEFT'
    layout = Table([[left_table, image]], colWidths=[240, 300])
    layout.setStyle(TableStyle([
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('GRID', (0, 0), (-1, -1), 1, colors.white)
    ]))
    elements.append(layout)
    elements.append(Spacer(1, 10))

def add_pagebreak(elements):
    elements.append(PageBreak())

def build_pdf_report(
    this_plate,config, kit, chemistry, sample_df, subpool, color_by, plate_map,
    align_df, formatted_counts, qc_200umi, qc_500umi, knee_plot_path,
    hmap_paths, 
    cb_settings, cb_results, cb_knee_path,
    violin1, violin2, violin1_filt, violin2_filt,
    stacked_main, stacked_mult,
    barcode_df,
    adata_obs, adata_obs_filt,combined_obs,
    main_tissue, mult_tissue, 
    this_sampletype, 
    qc_min_counts, 
    qc_min_genes, 
    qc_max_counts, 
    qc_pct_counts_mt, 
    qc_doublet_score
):
    elements = []
    
    ### SECTION 0: Title and overview
    add_title(elements, "UCI/Caltech IGVF Dataset Pipeline Report")
    overview_data = [
        ["", this_plate],
        ["Kit", kit],
        ["Chemistry", f"Parse {chemistry}"],
        ["Sample type", this_sampletype],
        ["Genotypes", "\n".join([g for g in sample_df['Genotype'].unique() if "/" not in g])],
        ["Tissues", "\n".join(sample_df['Tissue'].unique())],
        ["Subpool", subpool],
        ["Exome capture subpool", "Yes" if "EX" in subpool else "No"],
        [f"Expected {this_sampletype}", "100,000" if kit == "WT" else "1,000,000"]
    ]
    table = Table(overview_data)
    table.setStyle(TableStyle([
    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
    ('BACKGROUND', (0, 1), (-1, -1), colors.white),
    ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
    ('ALIGN', (0, 1), (0, -1), 'LEFT'),
    ('FONTNAME', (1, 1), (-1, -1), 'Helvetica')
    ]))

    add_layout(elements, table, plate_map)
    elements.append(Spacer(1, 20))
    
    section_num = 1
    table_num = 1
    
    ### OPTIONAL SECTION: Multiplexing
    if 'well_type' in sample_df.columns and any(sample_df['well_type'] == 'Multiplexed'):
        add_section_header(elements, f"{section_num}. Sample multiplexing")
        section_num += 1
        add_paragraph(elements, f"{main_tissue} was the main tisue on the plate while {mult_tissue} was multiplexed.")

        multiplexed_df = pd.read_csv('plots/multiplexed_genotype_key.csv')
        add_table(elements, multiplexed_df, title=f"Table {table_num}: Genotype multiplexing key")
        table_num += 1
        elements.append(Spacer(1, 10))
        

    ### SECTION 1: Alignment
    add_section_header(elements, f"{section_num}. Alignment")
    section_num += 1
    add_paragraph(elements, f"A total of {formatted_counts[2]} reads were aligned using GENCODE vM32 with mm39 assembly.")
    elements.append(Spacer(1, 10))
    add_table(elements, align_df, title=f"Table {table_num}: Pseudoalignment stats")
    table_num += 1
    
    ### SECTION 2: Cell recovery
    add_section_header(elements, f"{section_num}. Cell recovery")
    section_num += 1
    add_paragraph(elements, f"A total of {formatted_counts[0]} {this_sampletype} in {subpool} have >=200 UMIs and {formatted_counts[1]} {this_sampletype} have >= 500 UMIs.")
    add_image(elements, knee_plot_path, width = 375, height = 225)
    add_table(elements, qc_200umi, title= f"Table {table_num}: QC stats for {this_sampletype} >=200 UMI")
    table_num += 1
    add_table(elements, qc_500umi, title= f"Table {table_num}: QC stats for {this_sampletype} >=500 UMI")
    table_num += 1
    add_paragraph(elements, f"In this combinatorial barcoding assay, {this_sampletype} are barcoded in 96-well plates over three rounds. The first round assigns sample-specific barcodes, while the next two are applied to the pooled mix. Even distribution in the first round is important to ensure consistent sample representation.")
    add_image(elements, hmap_paths[0], width = 400, height = 220)
    add_image(elements, hmap_paths[1], width = 400, height = 220)
    add_image(elements, hmap_paths[2], width = 400, height = 220)
    
    ### SECTION 3: Background removal
    add_section_header(elements, f"{section_num}. Background removal")
    section_num += 1
    add_table(elements, cb_settings.drop_duplicates(), title= f"Table {table_num}: Cellbender parameters")
    table_num += 1
    add_table(elements, cb_results, title= f"Table {table_num}: Cellbender stats")
    table_num += 1
    add_image(elements, cb_knee_path, width = 425, height = 225)
    add_image(elements, violin1, width = 380, height = 150)
    add_image(elements, violin2, width = 380, height = 150)
    
    ### SECTION 4: Barcode sample map
    add_section_header(elements, f"{section_num}. Barcode sample map")
    section_num += 1
    add_paragraph(elements, f"The first round corresponds to sample barcoding. In Parse Bio assays, two barcodes are included in the first round: a random hexamer (R) and oligo-dT (T) barcode per well*. This reduces 3-prime bias and allows for the capture of full-length transcripts. These barcodes can be found in the barcode-1 region in the seqspec. For the {kit} kit, samples are distributed across {len(barcode_df['barcode'][barcode_df['parse barcode type'] == 'T'].unique().tolist())} wells in round 1. The following table includes our lab sample ID and IGVF sample accessions.")
    add_table(elements, barcode_df, title= f"Table {table_num}: Barcode sample map")
    table_num += 1
    add_pagebreak(elements)
    
    ### SECTION 5: Data processing workflow
    add_section_header(elements, f"{section_num}. Description of data processing workflow and h5ad files")
    section_num += 1
    add_paragraph(elements, f"Data from {subpool} were combined with all other subpools into one AnnData object per tissue as well as matching data from other plates for processing and cell type annotation. The AnnData contains detailed sample and cell metadata in the observation (obs) table and gene information in the variables table (var). Cell IDs (adata.obs.index) are re-formatted to be human-readable and unique across subpools and experiments by appending subpool and experiment IDs. The final AnnData contains both raw and CellBender denoised counts matrices in separate layers: adata.layers['raw_counts'] and adata.layers['cellbender_counts']. The current adata.X points to the CellBender denoised counts. Cells or nuclei are filtered by >0.5 cell_probability from CellBender analysis.")
    add_paragraph(elements, f"In addition to CellBender filtering, {this_sampletype} were also filtered by the following QC parameters:")
    qc_thresholds = pd.DataFrame({
        "Parameter": [
            "Minimum # UMIs",
            "Maximum # UMIs",
            "Minimum # genes",
            "Maximum % mito.",
            "Maximum doublet score"
        ],
        "Threshold": [
            qc_min_counts,
            qc_max_counts,
            qc_min_genes,
            qc_pct_counts_mt,
            qc_doublet_score
        ]
    })
    add_table(elements, qc_thresholds, title=f"Table {table_num}: QC thresholds")
    table_num += 1
    add_image(elements, violin1_filt, width = 380, height = 150)
    add_image(elements, violin2_filt, width = 380, height = 150)
    n_clust_main_tissue = combined_obs['leiden'][combined_obs['Tissue'] == main_tissue].astype(int).max() + 1
    n_clust_mult_tissue = combined_obs['leiden'][combined_obs['Tissue'] == mult_tissue].astype(int).max() + 1
    add_paragraph(elements, f"After filtering, data for each tissue was normalized and clustered using CellBender denoised counts. Normalization includes regression of technical parameters (number of genes expressed per cell and percent mitochondrial gene expression). Harmony was used to integrate exome capture and non-exome capture subpools for annotation. A Leiden clustering resolution of 1 was used for {n_clust_main_tissue} total clusters in {main_tissue} and {n_clust_mult_tissue} clusters in {mult_tissue}, using all tissue data across experiments. Cell type annotations were assigned per Leiden cluster or subcluster (leiden_R) using expression of canonical cell type marker genes.")
    add_image(elements, stacked_main, width = 450, height = 230)
    add_image(elements, stacked_mult, width = 450, height = 230)
    obs_desc = pd.read_csv("obs_keys_description.csv")
    add_table(elements, obs_desc, title=f"Table {table_num}: Description of obs keys in h5ad file")
    table_num += 1
    
    ### SECTION 6: Package versions
    add_section_header(elements, f"{section_num}. Package versions")
    section_num += 1
    package_versions_df = pd.DataFrame([
        ["kallisto", "0.50.1"],
        ["bustools", "0.43.2"],
        ["kb", "0.28.2"],
        ["Python", "3.9.0"],
        ["snakemake", "7.32.0"],
        ["cellbender", "0.3.2"],
        ["scanpy", "1.10.2"],
        ["scrublet", "0.2.3"],
        ["anndata", "0.10.8"],
        ["pandas", "2.2.2"],
        ["numpy", "1.26.4"],
        ["harmonypy", "0.0.9"]
    ], columns=["package", "version"])
    add_table(elements, package_versions_df, title=f"Table {table_num}: Package versions")
    table_num += 1

    ### SECTION 7: Contact info
    add_section_header(elements, f"{section_num}. Contact information and GitHub")
    section_num += 1
    add_paragraph(elements, "Ali Mortazavi, ali.mortazavi@uci.edu")
    add_paragraph(elements, "Elisabeth Rebboah, erebboah@uci.edu")

    pipeline_text = "Pipeline repo: "
    pipeline_link = '<a href="https://github.com/mortazavilab/parse_pipeline" color="blue">https://github.com/mortazavilab/parse_pipeline</a>'
    add_paragraph(elements, pipeline_text + pipeline_link)

    code_text = "Code used to make this report: "
    code_link = f'<a href="https://github.com/mortazavilab/8cube_manuscript/blob/main/plate_report/Subpool_report_{this_plate}.ipynb" color="blue">https://github.com/mortazavilab/8cube_manuscript/blob/main/plate_report/Subpool_report_{this_plate}.ipynb</a>'
    add_paragraph(elements, code_text + code_link)

    analysis_text = "Code used to process and annotate the data for main tissue: "
    analysis_link = f'<a href="https://github.com/mortazavilab/8cube_manuscript/blob/main/annotation/{main_tissue}.ipynb" color="blue">https://github.com/mortazavilab/8cube_manuscript/blob/main/annotation/{main_tissue}.ipynb</a>'
    add_paragraph(elements, analysis_text + analysis_link)

    analysis_text2 = "Code used to process and annotate the data for multiplexed tissue: "
    analysis_link2 = f'<a href="https://github.com/mortazavilab/8cube_manuscript/blob/main/annotation/{mult_tissue}.ipynb" color="blue">https://github.com/mortazavilab/8cube_manuscript/blob/main/annotation/{mult_tissue}.ipynb</a>'
    add_paragraph(elements, analysis_text2 + analysis_link2)

    
    return elements