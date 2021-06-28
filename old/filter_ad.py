import scrublet as scr
import scanpy as sc
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

aces = []
c1 = c2 = c3 = 0
for i in range(8):
    print(i)
    fname = "input_sce_Batch%d.h5ad" % (i + 1)
    adata = ad.read_h5ad(fname)
    c1 += sum(adata.obs.Phenotype != "BD")
    # Calculate QC metrics
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
    ]
    c2 += sum(adata.obs.Phenotype != "BD")
    # Remove doublets
    counts_matrix = adata.X
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs["scrublet_doublet_scores"] = doublet_scores
    adata.obs["scrublet_is_doublet"] = predicted_doublets
    adata = adata[
        adata.obs["scrublet_is_doublet"] == False,
    ]
    c3 += sum(adata.obs.Phenotype != "BD")
    aces.append(adata)
#    adata.write("input_sce_Batch%d_annotated.h5ad" % (i + 1))

print(c1)  # 560,020
print(c2)  # 462,229
print(c3)  # 454,622

merged_ace = ad.concat(aces)
merged_ace.obs_names_make_unique()  # (643988, 33538)


sc.pp.filter_genes(
    merged_ace,
    min_counts=None,
    min_cells=round(0.001 * merged_ace.shape[0]),
    max_counts=None,
    max_cells=None,
    inplace=True,
    copy=False,
)  # (643988, 24888)

merged_ace.X = csr_matrix(merged_ace.X)

"""
sc.pp.normalize_total(merged_ace)
sc.pp.log1p(merged_ace)

sc.pp.highly_variable_genes(merged_ace)
"""

merged_ace.write_h5ad("Ruzika_doublet_removed_all_counts.h5ad")

SZ_ace = merged_ace[
    merged_ace.obs.Phenotype != "BD",
]  # (454622, 24888)
SZ_ace.write_h5ad("Ruzika_doublet_removed_SZplusCON_counts.h5ad")

