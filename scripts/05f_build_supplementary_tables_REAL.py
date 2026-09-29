"""
Rebuilds Supplementary_Tables_REAL.xlsx from ONLY the real, verified files
produced in this analysis. Does NOT include Cell_Type_Proportions or
Microglia_Analysis sheets — those were fabricated in the old workbook and
were never actually computed for real. If you want them, they need to be
built from scratch and verified the same way the DE pipeline was.
"""
import pandas as pd

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'
OUT = f'{AD_PROJECT_DIR}/Supplementary_Tables_REAL.xlsx'

sheets = {
    'Sample_Metadata': pd.read_csv(f'{AD_PROJECT_DIR}/pseudobulk_sample_metadata.csv'),
    'All_Significant_DEGs': pd.read_csv(f'{AD_PROJECT_DIR}/deg_significant_padj05_REAL_annotated.csv'),
    'Cross_Dataset_Consistent_DEGs': pd.read_csv(f'{AD_PROJECT_DIR}/deg_high_confidence_cross_dataset_REAL.csv'),
    'Cross_Dataset_Consistency_Detail': pd.read_csv(f'{AD_PROJECT_DIR}/cross_dataset_consistency_REAL.csv'),
    'Pathway_Enrichment_Full': pd.read_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_full_REAL.csv'),
    'Pathway_Enrichment_Upregulated': pd.read_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_upregulated_REAL.csv'),
    'GSE296027_Standalone_CAVEAT': pd.read_csv(f'{AD_PROJECT_DIR}/deg_results_GSE296027_standalone.csv'),
}

# Try to add the downregulated enrichment sheet only if it has content
# (it may be empty — 0/49 terms significant — which is a real result, not
# a bug, so handle that gracefully rather than erroring).
try:
    down_df = pd.read_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_downregulated_REAL.csv')
    if len(down_df) > 0:
        sheets['Pathway_Enrichment_Downregulated'] = down_df
    else:
        print("Note: pathway_enrichment_downregulated_REAL.csv is empty (0 significant "
              "terms) — this is the correct, real result (underpowered 8-gene set), "
              "not an error. Adding a placeholder sheet noting this explicitly.")
        sheets['Pathway_Enrichment_Downregulated'] = pd.DataFrame(
            {'note': ['No terms significant after background correction (0/49 tested). '
                      'Downregulated gene set (n=8) too small for reliable enrichment.']})
except FileNotFoundError:
    print("pathway_enrichment_downregulated_REAL.csv not found — skipping that sheet.")

with pd.ExcelWriter(OUT, engine='openpyxl') as writer:
    for name, df in sheets.items():
        df.to_excel(writer, sheet_name=name[:31], index=False)  # Excel sheet name limit

# Add a README sheet explaining what's NOT here and why
readme = pd.DataFrame({
    'Note': [
        'This workbook replaces Supplementary_Tables.xlsx, which was found to be '
        'fabricated in its entirety.',
        '',
        'Sheets intentionally NOT included (were fabricated in the old workbook, '
        'never actually computed on real data):',
        '  - Cell_Type_Proportions',
        '  - Microglia_Analysis',
        '  - All_DEGs_Per_Dataset (per-cell-type — this analysis only ran whole-sample pseudobulk)',
        '  - KEGG_Enrichment / Reactome_Enrichment as separate sheets (real enrichment '
        'results are combined in Pathway_Enrichment_Full/Upregulated, with a Gene_set '
        'column distinguishing GO/KEGG/Reactome)',
        '',
        'GSE296027_Standalone_CAVEAT: this dataset was excluded from all directional '
        'DE/enrichment claims due to a complete confound between AD/Control status and '
        'sequencing plate. Its DE result is included here for transparency only — do '
        'not cite it as supporting evidence in the manuscript text.',
    ]
})
with pd.ExcelWriter(OUT, engine='openpyxl', mode='a') as writer:
    readme.to_excel(writer, sheet_name='README', index=False)

print(f"✓ Saved: {OUT}")
print(f"Sheets: {list(sheets.keys()) + ['README']}")
