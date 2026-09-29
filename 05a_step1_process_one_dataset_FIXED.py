"""
STEP 1 of 2 — process ONE dataset, save a small checkpoint file.
==================================================================

Run this once per dataset. Change GSE_TO_RUN below and re-run the cell
(or restart the Colab runtime first) for each of your datasets. Nothing
from a previous dataset needs to still be in memory — each run is
completely independent and writes its own small checkpoint file.

WHY THIS FIXES THE CRASH:
sc.read_h5ad() loads the *entire* clustered object: the counts you
need, but also the neighbor graph (.obsp, can be huge for 100k+ cells),
PCA/UMAP coordinates (.obsm), and cluster metadata (.uns) — none of
which pseudobulk aggregation uses. For GSE214244 (100k+ cells) that
extra weight is very likely what's blowing your RAM budget.

This script reads ONLY the three things it needs directly via h5py:
  - var index (gene names)
  - the counts layer (or .X if no counts layer), as sparse CSR
  - the sample_id obs column
It never constructs a full AnnData object, so obsp/obsm/uns are never
loaded at all.

OUTPUT (per dataset, small — a few MB, safe to write incrementally):
  /content/drive/MyDrive/AD_project/pseudobulk_parts/{GSE}_counts.csv
  /content/drive/MyDrive/AD_project/pseudobulk_parts/{GSE}_meta.csv

After you've run this once for every dataset you want included, run
05a_step2_combine_parts.py ONCE to merge everything into the final
pseudobulk_counts_matrix.csv / pseudobulk_sample_metadata.csv that the
R script expects. That combine step only reads small CSVs — no h5ad,
no scanpy, negligible RAM — so it's safe to (re)run any time, in a
fresh runtime, regardless of how many separate sessions you needed for
Step 1.
"""

import os
import re
import gc
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

warnings_shown = False

# =============================================================================
# >>> CHANGE THIS ONE LINE PER RUN <<<
# =============================================================================
GSE_TO_RUN = "GSE214244"
# =============================================================================

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'
PARTS_DIR = f'{AD_PROJECT_DIR}/pseudobulk_parts'
os.makedirs(PARTS_DIR, exist_ok=True)

# Crosswalks for the 2 datasets with no native sample_id column
gse143758_mapping = {
    'EZ.Batch1_HipR-AD-G1-4w': 'AD', 'EZ.Batch1_HipR-AD-G3-2w': 'AD',
    'EZ.Batch1_HipR-WT-G1-4w': 'Control', 'EZ.Batch1_HipR-WT-G3-2w': 'Control',
    'NP40.Batch1_Wt-Hip-S2-R': 'Control', 'NP40.Batch1_Untreated-Hip-S2-R': 'AD',
    'NP40.Batch3_Wt-Hip-S1-L': 'Control', 'NP40.Batch3_Wt-Hip-S2-L': 'Control',
    'NP40.Batch3_Untreated-Hip-S1-L': 'AD', 'NP40.Batch3_Untreated-Hip-S2-L': 'AD',
}
gse208683_mapping = {'WT_1': 'Control', 'WT_2': 'Control', 'AD_1': 'AD', 'AD_2': 'AD'}

DATASET_CONFIG = {
    'GSE140399': {'has_native_sample_id': True, 'crosswalk': None},
    'GSE140510': {'has_native_sample_id': True, 'crosswalk': None},
    'GSE143758': {'has_native_sample_id': False, 'crosswalk': gse143758_mapping},
    'GSE208683': {'has_native_sample_id': False, 'crosswalk': gse208683_mapping},
    'GSE214244': {'has_native_sample_id': True, 'crosswalk': None},
    'GSE296027': {'has_native_sample_id': True, 'crosswalk': None},
}

GSM_REGEX = re.compile(r'(GSM\d+)')


def decode_str_array(arr):
    """h5py returns bytes for string datasets in most anndata versions."""
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def read_h5_field(node):
    """Read an obs/var field (including the index), handling all three
    anndata on-disk string encodings seen in the wild:
      - plain h5py Dataset of bytes/str
      - categorical group: {'categories': [...], 'codes': [...]}
      - nullable/masked string array: {'values': [...], 'mask': [...]}
    Returns a numpy object array of python strings; masked/negative-code
    entries become None (caller should treat None as 'unmapped')."""
    if isinstance(node, h5py.Group):
        keys = set(node.keys())
        if {'categories', 'codes'} <= keys:
            cats = decode_str_array(node['categories'][:])
            codes = node['codes'][:]
            out = np.empty(len(codes), dtype=object)
            valid = codes >= 0
            out[valid] = cats[codes[valid]]
            out[~valid] = None
            return out
        elif {'values', 'mask'} <= keys:
            vals = node['values'][:]
            if vals.dtype.kind in ('S', 'O'):
                vals = decode_str_array(vals)
            vals = np.array(vals, dtype=object)
            mask = node['mask'][:].astype(bool)  # True = missing/null, per anndata convention
            frac_masked = mask.mean() if len(mask) else 0.0
            if frac_masked > 0.99:
                # A field where >99% of rows are "missing" is a red flag, not real
                # data — this is almost always a mask-convention/encoding mismatch
                # for THIS file, not 99% genuinely-null cells (an index in
                # particular can never legitimately be mostly-missing). Rather
                # than silently returning a column of Nones (which then fails
                # downstream as "0 matched" with no clue why), ignore the mask,
                # use the raw values, and say so loudly.
                print(f"  ⚠ WARNING: field mask flags {frac_masked:.1%} of rows as "
                      f"missing/null — that's almost certainly a mask-convention "
                      f"mismatch for this file, not real missingness. Ignoring the "
                      f"mask and using raw values instead. First 5 raw values: "
                      f"{list(vals[:5])}")
                return vals
            vals[mask] = None
            return vals
        else:
            raise ValueError(f"Unrecognized obs/var field encoding, group keys={sorted(keys)}")
    else:
        vals = node[:]
        if vals.dtype.kind in ('S', 'O'):
            vals = decode_str_array(vals)
        return vals


def read_obs_column(obs_grp, colname):
    return read_h5_field(obs_grp[colname])


def read_obs_index(obs_grp):
    idx_name = obs_grp.attrs.get('_index', '_index')
    if idx_name not in obs_grp:
        raise KeyError(f"Index field {idx_name!r} not found in .obs (keys: {list(obs_grp.keys())})")
    return read_h5_field(obs_grp[idx_name])


def read_gene_index(var_grp):
    idx_name = var_grp.attrs.get('_index', '_index')
    if idx_name not in var_grp:
        # very old anndata h5ad without an explicit _index attr
        idx_name = 'index' if 'index' in var_grp else list(var_grp.keys())[0]
    return decode_str_array(var_grp[idx_name][:])


def read_counts_matrix(f):
    """Return the counts layer (preferred) or .X as a scipy sparse CSR matrix,
    without ever materializing a dense copy or an AnnData object."""
    if 'layers' in f and 'counts' in f['layers']:
        grp = f['layers/counts']
        source = "layers['counts']"
    else:
        grp = f['X']
        source = ".X (WARNING: confirm this is raw counts, not log-normalized)"
    print(f"  Using {source}")

    if isinstance(grp, h5py.Group):  # sparse, stored as data/indices/indptr
        data = grp['data'][:]
        indices = grp['indices'][:]
        indptr = grp['indptr'][:]
        shape = grp.attrs.get('shape', None)
        enc = grp.attrs.get('encoding-type', 'csr_matrix')
        if shape is None:
            n_obs = f['obs'][list(f['obs'].keys())[0]].shape[0]
            n_var = f['var'][list(f['var'].keys())[0]].shape[0]
            shape = (n_obs, n_var)
        mat = sparse.csr_matrix((data, indices, indptr), shape=tuple(shape))
        if 'csc' in str(enc):
            mat = mat.tocsr()
    else:  # dense dataset
        mat = sparse.csr_matrix(grp[:])
    return mat


def process_dataset(gse):
    config = DATASET_CONFIG[gse]
    path = f'{AD_PROJECT_DIR}/{gse}/{gse}_clustered.h5ad'

    print(f"[{gse}] Opening {path} (h5py, minimal-RAM mode)...")
    if not os.path.exists(path):
        print(f"  ✗ File not found: {path}")
        return

    with h5py.File(path, 'r') as f:
        genes = read_gene_index(f['var'])
        print(f"  {len(genes)} genes in this dataset's panel")

        counts = read_counts_matrix(f)
        print(f"  Counts matrix: {counts.shape}")

        if config['has_native_sample_id']:
            if 'sample_id' not in f['obs']:
                print("  ✗ sample_id not found in .obs — skipping")
                return
            sample_id_raw = read_obs_column(f['obs'], 'sample_id').astype(str)

            master_meta_path = f'{AD_PROJECT_DIR}/master_metadata.csv'
            master_meta = pd.read_csv(master_meta_path)
            lookup = {(str(r['GSE']), str(r['GSM'])): r['condition_FINAL']
                      for _, r in master_meta.iterrows()}

            conditions = []
            for sid in sample_id_raw:
                m = GSM_REGEX.search(sid)  # search, not match — GSM may not be at position 0
                conditions.append(lookup.get((gse, m.group(1)), 'Unknown') if m else 'Unknown')
            conditions = np.array(conditions)
            group_ids = sample_id_raw
        else:
            barcodes = read_obs_index(f['obs'])
            crosswalk = config['crosswalk']
            group_ids = np.array([
                next((p for p in crosswalk if bc is not None and bc.startswith(p)), None)
                for bc in barcodes
            ], dtype=object)
            conditions = np.array([crosswalk.get(g, 'Unknown') if g else 'Unknown' for g in group_ids])

            n_matched = sum(g is not None for g in group_ids)
            if n_matched == 0:
                n_none_barcodes = sum(bc is None for bc in barcodes)
                print(f"  ✗ 0/{len(barcodes)} barcodes matched any crosswalk prefix.")
                print(f"    First 5 barcodes seen: {list(barcodes[:5])}")
                print(f"    Crosswalk prefixes expected: {list(crosswalk.keys())}")
                if n_none_barcodes == len(barcodes):
                    print(f"    -> ALL {len(barcodes)} barcodes came back as None from the "
                          f"index reader itself (this dataset's obs index is stored as a "
                          f"masked/nullable field). That's the actual problem — the crosswalk "
                          f"prefixes are probably fine. If you still see this after the mask-"
                          f"convention fix in read_h5_field, the index field is genuinely "
                          f"unreadable and needs manual inspection with h5py directly.")
                else:
                    print(f"    -> Barcodes ARE real strings, they just don't start with any "
                          f"crosswalk key. Fix the crosswalk keys above to match the real "
                          f"prefix format, then re-run. Not writing a checkpoint.")
                return

    valid_mask = np.isin(conditions, ['AD', 'Control'])
    n_valid = valid_mask.sum()
    print(f"  {n_valid} / {len(conditions)} cells resolved to AD/Control "
          f"({len(conditions) - n_valid} dropped as EXCLUDE_*/unmapped)")

    if n_valid == 0:
        print(f"  ✗ No usable cells for {gse} — check master_metadata.csv / crosswalk. Not writing a checkpoint.")
        return

    counts = counts[valid_mask]
    group_ids = group_ids[valid_mask]
    conditions = conditions[valid_mask]

    unique_samples = pd.unique(group_ids)
    print(f"  Aggregating to pseudobulk ({len(unique_samples)} samples)...")

    rows = {}
    meta_rows = []
    for sid in unique_samples:
        idx = (group_ids == sid)
        summed = np.asarray(counts[idx].sum(axis=0)).flatten().astype(np.int64)
        rows[f"{gse}_{sid}"] = pd.Series(summed, index=genes).groupby(level=0).sum()
        meta_rows.append({
            'sample': f"{gse}_{sid}", 'GSE': gse, 'sample_id': str(sid),
            'condition': conditions[idx][0], 'n_cells': int(idx.sum()),
            'total_UMI': int(summed.sum()),
        })

    part_counts_df = pd.DataFrame(rows)  # genes x samples-in-this-dataset
    part_meta_df = pd.DataFrame(meta_rows)

    out_counts = f'{PARTS_DIR}/{gse}_counts.csv'
    out_meta = f'{PARTS_DIR}/{gse}_meta.csv'
    part_counts_df.to_csv(out_counts)
    part_meta_df.to_csv(out_meta, index=False)
    print(f"  ✓ Saved checkpoint: {out_counts}  ({part_counts_df.shape})")
    print(f"  ✓ Saved checkpoint: {out_meta}")

    del counts, part_counts_df, part_meta_df, rows
    gc.collect()
    print(f"  ✓ Done with {gse}\n")


if __name__ == "__main__":
    if GSE_TO_RUN not in DATASET_CONFIG:
        raise ValueError(f"{GSE_TO_RUN} not in DATASET_CONFIG — add it first.")
    process_dataset(GSE_TO_RUN)
    print("Checkpoint files written for this dataset only.")
    print("Change GSE_TO_RUN above (restart runtime if you want a clean slate) "
          "and repeat for every remaining dataset.")
    print("Once ALL datasets you want are checkpointed in pseudobulk_parts/, "
          "run 05a_step2_combine_parts.py once — it's cheap, no scanpy/h5py needed.")
