#!/usr/bin/env python3
"""
LEAVE-ONE-DATASET-OUT VISUALIZATION
====================================
Creates publication-ready figures from the leave-one-out robustness analysis.

Input files (from 03_LEAVE_ONE_OUT_ROBUSTNESS.R):
  - loocv_results_summary.csv
  - loocv_consistency_table.csv
  - loocv_top_genes_persistence.csv

Output figures:
  - loocv_bar_plot.pdf (n_significant genes across runs)
  - loocv_consistency_heatmap.pdf (top genes persistence across runs)
  - loocv_summary_figure.pdf (comprehensive summary)

Usage:
    python 04_LOOCV_VISUALIZATION.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

OUTPUT_DIR = Path(".")

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_significant_genes_bar():
    """
    Bar plot: Number of significant genes across all leave-one-out runs.
    Shows that the number of significant genes is stable despite dropping datasets.
    """
    
    # Load data
    summary = pd.read_csv("loocv_results_summary.csv")
    
    # Add reference run (full analysis)
    # Read from the original full analysis if available
    reference = pd.DataFrame({
        'run': ['Full_5_datasets'],
        'dropped_dataset': ['None'],
        'n_genes_tested': [20629],  # From manuscript
        'n_significant': [301],      # From manuscript
        'n_upregulated': [290],
        'n_downregulated': [11],
        'top10_genes': ['Ccl3|Itgax|Ccl4|Lilrb4a|Cst7|Plaur|Gm26714|Ctse|Asb10|Tnfsf8']
    })
    
    summary = pd.concat([reference, summary], ignore_index=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot counts
    x_pos = np.arange(len(summary))
    width = 0.35
    
    ax.bar(x_pos - width/2, summary['n_upregulated'], width, 
           label='Upregulated', color='#e74c3c', alpha=0.8)
    ax.bar(x_pos + width/2, summary['n_downregulated'], width,
           label='Downregulated', color='#3498db', alpha=0.8)
    
    # Labels
    run_labels = [d.replace('Drop_', '').replace('_5_datasets', ' (Full)') 
                  for d in summary['dropped_dataset']]
    ax.set_xticks(x_pos)
    ax.set_xticklabels(run_labels, rotation=45, ha='right')
    
    ax.set_ylabel('Number of Significant Genes (FDR < 0.05)', fontsize=12)
    ax.set_xlabel('Analysis Run', fontsize=12)
    ax.set_title('Leave-One-Dataset-Out Robustness: Significant Gene Counts', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'loocv_bar_plot.pdf', dpi=300, bbox_inches='tight')
    print("✓ Saved: loocv_bar_plot.pdf")
    plt.close()

def plot_top_genes_persistence():
    """
    Heatmap: Showing whether top 10 genes from full analysis appear in each leave-one-out run.
    """
    
    # Load data
    persistence = pd.read_csv("loocv_top_genes_persistence.csv")
    
    # Extract presence/absence from text like "Yes (rank 1)"
    heatmap_data = pd.DataFrame(index=persistence['gene'])
    
    for col in persistence.columns[1:]:  # Skip gene column
        heatmap_data[col.replace('Drop_', '')] = persistence[col].apply(
            lambda x: 1 if 'Yes' in str(x) else 0
        )
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 6))
    
    sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='RdYlGn', 
                cbar_kws={'label': 'Present in Run'}, ax=ax, linewidths=0.5)
    
    ax.set_title('Top 10 Genes: Presence Across Leave-One-Dataset-Out Runs', 
                 fontsize=14, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=12)
    ax.set_xlabel('Excluded Dataset', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'loocv_consistency_heatmap.pdf', dpi=300, bbox_inches='tight')
    print("✓ Saved: loocv_consistency_heatmap.pdf")
    plt.close()

def plot_summary_stats():
    """
    Summary figure: Multiple panels showing robustness metrics.
    """
    
    # Load data
    consistency = pd.read_csv("loocv_consistency_table.csv")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Panel A: Significant genes in each run
    ax = axes[0]
    dropped = ['Dropped: ' + d for d in consistency['dropped_dataset']]
    ax.barh(dropped, consistency['n_sig_in_run'], color='#2ecc71', alpha=0.8)
    ax.set_xlabel('Number of Significant Genes', fontsize=11)
    ax.set_title('A) Significant Genes per Run\n(FDR < 0.05)', fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    # Add reference line for full analysis (301 genes)
    ax.axvline(301, color='red', linestyle='--', linewidth=2, label='Full analysis (301 genes)')
    ax.legend()
    
    # Panel B: Retention of top genes
    ax = axes[1]
    ax.bar(dropped, consistency['pct_top_genes_retained'], color='#3498db', alpha=0.8)
    ax.set_ylabel('% of Top 10 Genes Retained', fontsize=11)
    ax.set_title('B) Persistence of Top 10 Genes\n(from full analysis)', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 110])
    ax.grid(axis='y', alpha=0.3)
    
    # Add reference line at 100%
    ax.axhline(100, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Full retention')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'loocv_summary_figure.pdf', dpi=300, bbox_inches='tight')
    print("✓ Saved: loocv_summary_figure.pdf")
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("\n" + "="*70)
    print("LEAVE-ONE-OUT ROBUSTNESS VISUALIZATION")
    print("="*70)
    
    # Check if input files exist
    required_files = [
        "loocv_results_summary.csv",
        "loocv_consistency_table.csv",
        "loocv_top_genes_persistence.csv"
    ]
    
    missing = [f for f in required_files if not Path(f).exists()]
    if missing:
        print(f"\n✗ Missing input files: {', '.join(missing)}")
        print(f"  (Run 03_LEAVE_ONE_OUT_ROBUSTNESS.R first)")
        return
    
    print("\n✓ All input files found")
    
    # Create figures
    print("\nCreating figures...")
    try:
        plot_significant_genes_bar()
        plot_top_genes_persistence()
        plot_summary_stats()
    except Exception as e:
        print(f"✗ Error creating figures: {e}")
        return
    
    print("\n" + "="*70)
    print("✅ VISUALIZATION COMPLETE")
    print("="*70)
    print("\nGenerated figures:")
    print("  - loocv_bar_plot.pdf")
    print("  - loocv_consistency_heatmap.pdf")
    print("  - loocv_summary_figure.pdf")
    print("\nFor manuscript:")
    print("  - Include loocv_summary_figure.pdf as supplementary figure")
    print("  - Reference in Results section: 'Leave-one-dataset-out analysis'")

if __name__ == "__main__":
    main()
