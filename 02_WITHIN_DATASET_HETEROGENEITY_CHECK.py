#!/usr/bin/env python3
"""
WITHIN-DATASET HETEROGENEITY CHECK
===================================
Examines master_metadata.csv to identify:
- GSE214244: EC (entorhinal cortex) vs DCN (deep cerebellar nuclei) samples
- GSE143758: EZ vs NP40 protocols, 2-week vs 4-week timepoints
- Balance across conditions (AD vs Control)

This helps identify potential confounders not controlled in the pseudobulk DE model.

Usage:
    python 02_WITHIN_DATASET_HETEROGENEITY_CHECK.py

Outputs:
    - heterogeneity_report.txt (for Limitations section)
    - heterogeneity_details.csv (detailed breakdown)
"""

import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

METADATA_FILE = Path("/path/to/master_metadata.csv")  # UPDATE THIS

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def check_gse214244(df):
    """
    Check GSE214244 for EC (entorhinal cortex) vs DCN (deep cerebellar nuclei) split.
    """
    print("\n" + "="*70)
    print("GSE214244: Region/Protocol Heterogeneity Check")
    print("="*70)
    
    subset = df[df['GSE'] == 'GSE214244'].copy()
    
    if len(subset) == 0:
        print("⚠ No GSE214244 samples found in metadata")
        return None
    
    print(f"\nTotal samples: {len(subset)}")
    
    # Look for region/protocol indicators in sample names or columns
    print("\nSample names (first 10):")
    for i, name in enumerate(subset['sample_id'].head(10)):
        print(f"  {i+1}. {name}")
    
    # Check if brain_region_FINAL column exists
    if 'brain_region_FINAL' in subset.columns:
        regions = subset['brain_region_FINAL'].value_counts()
        print(f"\nRegions in data:")
        for region, count in regions.items():
            print(f"  {region}: {count} samples")
        
        # Check balance by condition
        print(f"\nRegion × Condition crosstab:")
        crosstab = pd.crosstab(subset['brain_region_FINAL'], subset['condition_FINAL'], margins=True)
        print(crosstab)
        
        # Identify potential confounding
        ec_samples = subset[subset['brain_region_FINAL'].str.contains('EC|entorhinal', case=False, na=False)]
        dcn_samples = subset[subset['brain_region_FINAL'].str.contains('DCN|cerebel', case=False, na=False)]
        
        if len(ec_samples) > 0 and len(dcn_samples) > 0:
            print(f"\n⚠️  HETEROGENEITY DETECTED:")
            print(f"  - EC samples: {len(ec_samples)} (AD: {(ec_samples['condition_FINAL']=='AD').sum()}, Control: {(ec_samples['condition_FINAL']=='Control').sum()})")
            print(f"  - DCN samples: {len(dcn_samples)} (AD: {(dcn_samples['condition_FINAL']=='AD').sum()}, Control: {(dcn_samples['condition_FINAL']=='Control').sum()})")
            
            if (ec_samples['condition_FINAL'] == 'AD').sum() == len(ec_samples) or \
               (dcn_samples['condition_FINAL'] == 'AD').sum() == len(dcn_samples):
                print(f"\n  🚨 POTENTIAL CONFOUND: One region is exclusively one condition!")
                return "CONFOUNDED"
            else:
                print(f"\n  ℹ  Regions are mixed by condition, but should be noted as a source of variation.")
                return "MIXED"
    
    else:
        print("\n⚠ 'brain_region_FINAL' column not found in metadata")
        print("  Checking 'notes' or 'sample_name' for EC/DCN indicators...")
        
        ec_count = subset['sample_id'].str.contains('EC|entorhinal', case=False, na=False).sum()
        dcn_count = subset['sample_id'].str.contains('DCN|cerebel', case=False, na=False).sum()
        
        if ec_count > 0 or dcn_count > 0:
            print(f"  EC-like samples: {ec_count}")
            print(f"  DCN-like samples: {dcn_count}")
            return "DETECTED_IN_NAMES"
    
    return None

def check_gse143758(df):
    """
    Check GSE143758 for EZ vs NP40 protocols and 2w vs 4w timepoints.
    """
    print("\n" + "="*70)
    print("GSE143758: Protocol and Timepoint Heterogeneity Check")
    print("="*70)
    
    subset = df[df['GSE'] == 'GSE143758'].copy()
    
    if len(subset) == 0:
        print("⚠ No GSE143758 samples found in metadata")
        return None
    
    print(f"\nTotal samples in metadata: {len(subset)}")
    print(f"Note: Only 10/37 samples are in the loaded hippocampus subset")
    
    # Look for protocol/timepoint indicators
    print("\nSample names (first 15):")
    for i, name in enumerate(subset['sample_id'].head(15)):
        print(f"  {i+1}. {name}")
    
    # Check for EZ vs NP40
    ez_count = subset['sample_id'].str.contains('EZ|ez', case=False, na=False).sum()
    np40_count = subset['sample_id'].str.contains('NP40|np40', case=False, na=False).sum()
    
    timepoint_2w = subset['sample_id'].str.contains('2w|2W|2-week', case=False, na=False).sum()
    timepoint_4w = subset['sample_id'].str.contains('4w|4W|4-week', case=False, na=False).sum()
    
    print(f"\nProtocols detected in sample names:")
    print(f"  EZ (TRIzol): {ez_count} samples")
    print(f"  NP40 (lysis buffer): {np40_count} samples")
    
    if ez_count > 0 and np40_count > 0:
        print(f"\n⚠️  PROTOCOL HETEROGENEITY DETECTED:")
        print(f"  - Dataset mixes {ez_count} EZ and {np40_count} NP40 samples")
        print(f"  - These protocols can produce different biases (e.g., RNA degradation, cell-type recovery)")
        
        # Check if confounded with condition
        ez_samples = subset[subset['sample_id'].str.contains('EZ|ez', case=False, na=False)]
        np40_samples = subset[subset['sample_id'].str.contains('NP40|np40', case=False, na=False)]
        
        print(f"\n  EZ × Condition:")
        print(f"    - {(ez_samples['condition_FINAL']=='AD').sum()} AD, {(ez_samples['condition_FINAL']=='Control').sum()} Control")
        print(f"\n  NP40 × Condition:")
        print(f"    - {(np40_samples['condition_FINAL']=='AD').sum()} AD, {(np40_samples['condition_FINAL']=='Control').sum()} Control")
        
        heterogeneity_type = "PROTOCOL_MIXED"
    else:
        print(f"\n✓ Only one protocol detected (or unclear from names)")
        heterogeneity_type = None
    
    # Check timepoints
    if timepoint_2w > 0 or timepoint_4w > 0:
        print(f"\nTimepoints detected in sample names:")
        print(f"  2-week: {timepoint_2w} samples")
        print(f"  4-week: {timepoint_4w} samples")
        
        if timepoint_2w > 0 and timepoint_4w > 0:
            print(f"\n⚠️  TIMEPOINT HETEROGENEITY DETECTED:")
            print(f"  - Dataset mixes {timepoint_2w} 2-week and {timepoint_4w} 4-week samples")
            
            w2_samples = subset[subset['sample_id'].str.contains('2w|2W|2-week', case=False, na=False)]
            w4_samples = subset[subset['sample_id'].str.contains('4w|4W|4-week', case=False, na=False)]
            
            print(f"\n  2-week × Condition:")
            print(f"    - {(w2_samples['condition_FINAL']=='AD').sum()} AD, {(w2_samples['condition_FINAL']=='Control').sum()} Control")
            print(f"\n  4-week × Condition:")
            print(f"    - {(w4_samples['condition_FINAL']=='AD').sum()} AD, {(w4_samples['condition_FINAL']=='Control').sum()} Control")
            
            if heterogeneity_type:
                heterogeneity_type += "_AND_TIMEPOINT"
            else:
                heterogeneity_type = "TIMEPOINT_MIXED"
    
    return heterogeneity_type

def check_all_datasets(df):
    """
    Scan all datasets for unbalanced covariates that might confound the results.
    """
    print("\n" + "="*70)
    print("GENERAL COVARIATE BALANCE CHECK (All Datasets)")
    print("="*70)
    
    for gse in df['GSE'].unique():
        subset = df[df['GSE'] == gse]
        usable = subset[subset['usable_for_DE'] == 'yes']
        
        if len(usable) > 0:
            print(f"\n{gse}:")
            print(f"  Total samples in data: {len(subset)}")
            print(f"  Usable for DE: {len(usable)}")
            
            # Check condition balance
            cond_counts = usable['condition_FINAL'].value_counts()
            print(f"  Conditions: {dict(cond_counts)}")
            
            # Check for uneven sample distribution in brain regions (if available)
            if 'brain_region_FINAL' in usable.columns:
                region_cond = pd.crosstab(usable['brain_region_FINAL'], usable['condition_FINAL'])
                if len(region_cond) > 1:
                    print(f"  Region imbalance:")
                    for region in region_cond.index:
                        print(f"    {region}: {dict(region_cond.loc[region])}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("\n" + "="*70)
    print("WITHIN-DATASET HETEROGENEITY CHECK")
    print("="*70)
    print(f"\nMetadata file: {METADATA_FILE}")
    
    # Load metadata
    try:
        metadata = pd.read_csv(METADATA_FILE)
        print(f"✓ Loaded metadata: {len(metadata)} rows")
        print(f"  Columns: {list(metadata.columns)[:5]}... (showing first 5)")
    except Exception as e:
        print(f"✗ Could not load metadata: {e}")
        return None
    
    # Run checks
    result_214244 = check_gse214244(metadata)
    result_143758 = check_gse143758(metadata)
    
    check_all_datasets(metadata)
    
    # Generate report for Limitations
    print("\n" + "="*70)
    print("LIMITATIONS SECTION ADDITIONS")
    print("="*70)
    
    report = "\n## Within-Dataset Heterogeneity\n\n"
    
    if result_214244:
        if result_214244 == "CONFOUNDED":
            report += "**GSE214244** includes samples from two anatomically distinct regions (entorhinal cortex and deep cerebellar nuclei) that are not balanced across conditions. This confounding was not explicitly controlled in the pseudobulk DE model and may inflate region-specific effects attributed to AD status.\n\n"
        else:
            report += "**GSE214244** includes samples from two anatomically distinct regions (entorhinal cortex and deep cerebellar nuclei). While these are mixed by condition, regional transcriptional differences independent of AD status may contribute to the pooled pseudobulk signal.\n\n"
    
    if result_143758:
        if "PROTOCOL" in result_143758:
            report += "**GSE143758** used two different sample preparation protocols (EZ/TRIzol and NP40 lysis). These protocols can produce different biases in RNA recovery and cell-type composition; the protocol effect was not explicitly modeled as a batch covariate in the pseudobulk DE design.\n\n"
        if "TIMEPOINT" in result_143758:
            report += "**GSE143758** combines samples from two timepoints (2-week and 4-week). Timepoint-specific transcriptional changes may contribute to the overall DE signal independently of AD status.\n\n"
    
    report += """Additionally, only 10 of 37 curated samples for GSE143758 are represented in the loaded hippocampus-restricted subset; results for this dataset should be interpreted as specific to this anatomical/temporal subset rather than the full original experiment."""
    
    # Save report
    with open('heterogeneity_report.txt', 'w') as f:
        f.write(report)
    
    print("\n" + report)
    print("\n✓ Saved to: heterogeneity_report.txt")
    
    # Save detailed results
    details = pd.DataFrame({
        'Dataset': ['GSE214244', 'GSE143758'],
        'Heterogeneity_Type': [result_214244 or "None detected", result_143758 or "None detected"],
        'Recommendation': [
            "Note in Methods 2.1 and Limitations; consider region × condition interaction if significant",
            "Note in Methods 2.1 and Limitations; protocol and timepoint effects may contribute to signal"
        ]
    })
    
    details.to_csv('heterogeneity_details.csv', index=False)
    print("✓ Saved to: heterogeneity_details.csv")
    
    return report

if __name__ == "__main__":
    report = main()
    print("\n" + "="*70)
    print("✅ HETEROGENEITY CHECK COMPLETE")
    print("="*70)
    print("\nNext steps:")
    print("  1. Review heterogeneity_report.txt")
    print("  2. Add relevant sections to Methods 2.1 (dataset description)")
    print("  3. Add to Limitations section")
