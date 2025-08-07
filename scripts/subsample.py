import scanpy as sc
import numpy as np
import sys

in_h5ad = sys.argv[1]
out_h5ad = in_h5ad.replace('.h5ad', '_subsampled20k.h5ad')

adata = sc.read_h5ad(in_h5ad)
print(f"Loaded: {in_h5ad}, {adata.shape}", flush=True)

if adata.n_obs >= 20000:
    np.random.seed(42)
    subsample_indices = np.random.choice(adata.n_obs, size=20000, replace=False)
    adata_sub = adata[subsample_indices].copy()
    adata_sub.write(out_h5ad)
    print(f"Saved: {out_h5ad}, {adata_sub.shape}", flush=True)
else:
    print(f"Skipping: only {adata.n_obs} cells", flush=True)

