#!/usr/bin/env python3
"""make_srvf_data.py -- SRVF (elastic) tangent-PCA characters for the tree.

Parallel to make_revbayes_data.R, but the continuous characters are SRVF
tangent-PCA scores instead of coordinate-PCA. SRVF is the faithful elastic-shape
representation: on these 39 pots its inter-type distances recover the true SRVF
Karcher distance matrix at Spearman 0.98 (coordinate-PCA only 0.73).

Uses fdasrsf (same package/settings as the corpus Karcher pipeline:
mode='O' open curves, no scaling, no rotation).

Run:  python3 make_srvf_data.py 20        # 20 SRVF PCs (default)
Writes data/mini_character_srvf_<N>.nex (fixed-point, INTERLEAVE=NO).
Reuses data/mini_taxa.tsv (ages are identical -- only the characters change).
"""
import sys, os, numpy as np, pandas as pd
from fdasrsf.curve_stats import fdacurve

N_PC = int(sys.argv[1]) if len(sys.argv) > 1 else 20
HERE = "/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
DATA = "/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"

sel  = pd.read_csv(f"{HERE}/selected_taxa.csv")
prof = pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv")
types = list(sel.type)

# beta: (n_dim=2, T=250, n_curves) in fixed point order
beta = np.stack([prof[prof.type == t].sort_values("point_order")[["x", "y"]]
                 .to_numpy().T for t in types], axis=2)

obj = fdacurve(beta, mode='O', N=beta.shape[1], scale=False)
obj.karcher_mean(rotation=False)
obj.srvf_align(rotation=False)
obj.shape_pca(no=max(N_PC, 20))
S = np.asarray(obj.coef)[:, :N_PC]                       # (39, N_PC) tangent-PCA scores
ve = np.asarray(obj.s); ve = ve / ve.sum()
print(f"SRVF tangent-PCA: {N_PC} PCs, {100*np.cumsum(ve)[N_PC-1]:.2f}% cumulative variance "
      f"(PC1 alone {100*ve[0]:.1f}%)")

# --- write NEXUS (pilot format, fixed-point; RevBayes rejects scientific notation)
os.makedirs(f"{HERE}/data", exist_ok=True)
nex = f"{HERE}/data/mini_character_srvf_{N_PC}.nex"
pad = max(len(t) for t in types) + 4
with open(nex, "w") as f:
    f.write("#NEXUS\n[SRVF tangent-PCA, make_srvf_data.py]\nBEGIN DATA;\n")
    f.write(f"  DIMENSIONS NTAX={len(types)} NCHAR={N_PC};\n")
    f.write("  FORMAT DATATYPE=CONTINUOUS MISSING=? GAP=- INTERLEAVE=NO;\n  MATRIX\n")
    for t, row in zip(types, S):
        f.write("    " + t.ljust(pad) + "\t".join(f"{v:.12f}" for v in row) + "\n")
    f.write("  ;\nEND;\n")

# sanity: taxa match the shared taxa file
tsv = pd.read_csv(f"{HERE}/data/mini_taxa.tsv", sep="\t")
assert set(types) == set(tsv.taxon), "taxon set mismatch vs mini_taxa.tsv"
assert not np.any([ "e" in f"{v:.12f}" for v in S.flatten() ]), "scientific notation leaked"
print(f"wrote {nex}  (NTAX={len(types)} NCHAR={N_PC}); taxa match mini_taxa.tsv")
