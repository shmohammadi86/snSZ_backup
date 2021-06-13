import scrublet as scr
import scanpy as sc
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

adata = ad.read_h5ad("Panos_raw_sce_simplified.h5ad")  # (187277, 33226)
adata.var["mt"] = adata.var_names.str.startswith(
    "MT-"
)  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True,
)
adata = adata[
    (1000 <= adata.obs.n_genes_by_counts)
    & (adata.obs.pct_counts_mt <= 10)
    & (adata.obs.total_counts <= 50000),
]  # (144543, 33226)

aces = []
sets = np.unique(adata.obs["set_ID"])
for set_id in sets:
    print(set_id)
    sub_adata = adata[
        adata.obs["set_ID"] == set_id,
    ]
    # Remove doublets
    counts_matrix = sub_adata.X
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    sub_adata.obs["scrublet_doublet_scores"] = doublet_scores
    sub_adata.obs["scrublet_is_doublet"] = predicted_doublets
    sub_adata = sub_adata[
        sub_adata.obs["scrublet_is_doublet"] == False,
    ]
    aces.append(sub_adata)

merged_ace = ad.concat(aces)  # (140917, 33226)
merged_ace.obs_names_make_unique()

sc.pp.filter_genes(
    merged_ace,
    min_counts=None,
    min_cells=round(0.001 * merged_ace.shape[0]),
    max_counts=None,
    max_cells=None,
    inplace=True,
    copy=False,
)  # (140917, 24458)

merged_ace.X = csr_matrix(merged_ace.X)

"""
sc.pp.normalize_total(merged_ace)
sc.pp.log1p(merged_ace)

sc.pp.highly_variable_genes(merged_ace)
"""

merged_ace.write_h5ad("Panos_doublet_removed_counts.h5ad")
