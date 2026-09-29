"""
Adds a confidence tier to cell_annotation_master_6datasets_CORRECTED.csv
and corrects the 4 clusters that got a specific cell-type label despite
having ZERO matched marker genes -- those get demoted to
'Ambiguous_needs_manual_review' rather than silently kept as confident
calls. Does not touch clusters with real (even if weak) evidence.
"""
import pandas as pd

df = pd.read_csv('cell_annotation_master_6datasets_CORRECTED.csv')

def confidence_tier(row):
    n = row['new_match_count'] if pd.notna(row['new_match_count']) and row['new_match_count'] > 0 else row['match_count']
    if pd.isna(n) or n == 0:
        return 'NO_EVIDENCE'
    elif n == 1:
        return 'WEAK_1_gene'
    elif n <= 3:
        return 'MODERATE_2to3_genes'
    else:
        return 'STRONG_4plus_genes'

df['annotation_confidence'] = df.apply(confidence_tier, axis=1)

# Demote the 4 zero-evidence clusters that got a SPECIFIC label instead of
# falling back to Ambiguous/Unclassified like the other 53 zero-evidence
# clusters did.
zero_evidence_specific = (
    (df['annotation_confidence'] == 'NO_EVIDENCE') &
    (~df['final_cell_type'].isin(['Neurons_Unclassified', 'Ambiguous_needs_manual_review']))
)
print(f"Demoting {zero_evidence_specific.sum()} clusters with a specific label but zero evidence:")
print(df.loc[zero_evidence_specific, ['GSE', 'cluster', 'n_cells', 'final_cell_type']].to_string())

df['final_cell_type_ORIGINAL'] = df['final_cell_type']  # keep for transparency
df.loc[zero_evidence_specific, 'final_cell_type'] = 'Ambiguous_needs_manual_review'

df.to_csv('cell_annotation_master_6datasets_CONFIDENCE_TIERED.csv', index=False)
print(f"\n✓ Saved: cell_annotation_master_6datasets_CONFIDENCE_TIERED.csv")
print(f"\nConfidence tier distribution:")
print(df['annotation_confidence'].value_counts())
print(f"\nFor the manuscript: report cell-type proportions using only STRONG_4plus_genes and")
print(f"MODERATE_2to3_genes tiers as primary evidence; note WEAK_1_gene and NO_EVIDENCE/Ambiguous")
print(f"clusters explicitly as lower-confidence in the methods or a supplementary caveat.")
