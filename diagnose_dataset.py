"""
DIAGNOSTIC — run this before patching sample-matching logic again.

Prints, for one dataset:
  - every column available in .obs (so we can see if there's a
    sample/batch/orig.ident-style column we should be using instead of
    guessing at barcode prefixes)
  - the first 20 raw obs index values (barcodes) so we can see their
    actual format
  - the first 20 values of every non-numeric-looking obs column, so we
    can see if any of them already contain something crosswalk-able
    (e.g. 'WT_1', 'AD_1', or a GSM id)

No counts are loaded. This is metadata-only and should run in seconds
regardless of dataset size.
"""

import h5py
import numpy as np

GSE_TO_INSPECT = "GSE208683"   # change as needed
AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'
path = f'{AD_PROJECT_DIR}/{GSE_TO_INSPECT}/{GSE_TO_INSPECT}_clustered.h5ad'


def decode_str_array(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


with h5py.File(path, 'r') as f:
    print(f"=== {GSE_TO_INSPECT} ===")
    print(f"Top-level keys: {list(f.keys())}")

    obs = f['obs']
    print(f"\n.obs attrs: {dict(obs.attrs)}")
    print(f".obs keys: {list(obs.keys())}")

    idx_name = obs.attrs.get('_index', '_index')
    print(f"\nIndex column name: {idx_name!r}")
    if idx_name in obs:
        idx_node = obs[idx_name]
        if isinstance(idx_node, h5py.Group):
            cats = decode_str_array(idx_node['categories'][:])
            codes = idx_node['codes'][:20]
            vals = cats[codes]
        else:
            vals = idx_node[:20]
            if vals.dtype.kind in ('S', 'O'):
                vals = decode_str_array(vals)
        print(f"First 20 index values (barcodes):\n{list(vals)}")
    else:
        print(f"  ✗ {idx_name!r} not found directly in .obs — check keys above")

    print("\n--- other obs columns (first 10 values each) ---")
    for col in obs.keys():
        if col == idx_name:
            continue
        node = obs[col]
        try:
            if isinstance(node, h5py.Group):
                if 'categories' in node:
                    cats = decode_str_array(node['categories'][:])
                    codes = node['codes'][:10]
                    vals = cats[codes]
                    print(f"  [{col}] (categorical, {len(cats)} categories): {list(vals)}")
                else:
                    print(f"  [{col}] is a group with keys {list(node.keys())} — skipping preview")
            else:
                vals = node[:10]
                if vals.dtype.kind in ('S', 'O'):
                    vals = decode_str_array(vals)
                print(f"  [{col}] dtype={node.dtype}: {list(vals)}")
        except Exception as e:
            print(f"  [{col}] — could not preview ({e})")

    print("\n--- var (genes) ---")
    var = f['var']
    var_idx_name = var.attrs.get('_index', '_index')
    print(f"var index name: {var_idx_name!r}, n_genes: {var[var_idx_name].shape[0] if var_idx_name in var else 'unknown'}")
