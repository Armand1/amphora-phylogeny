#!/usr/bin/env python3
"""make_srvf_data_77taxon.py -- SRVF tangent-PCA characters + ages for the neckbelly set.
Self-contained sibling of make_srvf_data.py (reads selected_taxa_neckbelly.csv, writes both
the NEXUS and the taxa TSV). Same fdasrsf settings (mode='O', no scale/rotation).
Run: python3 make_srvf_data_77taxon.py 38
"""
import sys, os, numpy as np, pandas as pd
from fdasrsf.curve_stats import fdacurve

N_PC = int(sys.argv[1]) if len(sys.argv) > 1 else 38
HERE = "/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
DATA = "/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"

sel  = pd.read_csv(f"{HERE}/selected_taxa_neckbelly.csv")
prof = pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv")
types = list(sel.type)

beta = np.stack([prof[prof.type == t].sort_values("point_order")[["x","y"]].to_numpy().T
                 for t in types], axis=2)
obj = fdacurve(beta, mode='O', N=beta.shape[1], scale=False)
obj.karcher_mean(rotation=False); obj.srvf_align(rotation=False)
obj.shape_pca(no=N_PC)
S = np.asarray(obj.coef).T[:, :N_PC]                      # (n_taxa, N_PC)
assert S.shape[0] == len(types), f"{S.shape[0]} != {len(types)}"
ve = np.asarray(obj.s); ve = ve/ve.sum()
print(f"SRVF tPCA on {len(types)} curves: {N_PC} PCs, PC1 {100*ve[0]:.1f}%, "
      f"cum {100*np.cumsum(ve)[N_PC-1]:.1f}%")

os.makedirs(f"{HERE}/data", exist_ok=True)
nex = f"{HERE}/data/mini_character_neckbelly_srvf_{N_PC}.nex"
pad = max(len(t) for t in types) + 4
with open(nex, "w") as f:
    f.write("#NEXUS\n[SRVF tPCA neckbelly, make_srvf_data_77taxon.py]\nBEGIN DATA;\n")
    f.write(f"  DIMENSIONS NTAX={len(types)} NCHAR={N_PC};\n")
    f.write("  FORMAT DATATYPE=CONTINUOUS MISSING=? GAP=- INTERLEAVE=NO;\n  MATRIX\n")
    for t, row in zip(types, S):
        f.write("    " + t.ljust(pad) + "\t".join(f"{v:.12f}" for v in row) + "\n")
    f.write("  ;\nEND;\n")
assert not any("e" in f"{v:.12f}" for v in S.flatten()), "sci notation leaked"

# taxa TSV (ages BP, datum 1950)
ages = prof[prof.type.isin(types)].drop_duplicates("type")[["type","date_earliest","date_latest"]]
ages = ages.assign(min_age=1950-ages.date_latest, max_age=1950-ages.date_earliest)
ages[["type","min_age","max_age"]].rename(columns={"type":"taxon"}).to_csv(
    f"{HERE}/data/mini_taxa_neckbelly.tsv", sep="\t", index=False)
print(f"wrote {os.path.basename(nex)} (NTAX={len(types)} NCHAR={N_PC}) + mini_taxa_neckbelly.tsv")
print(f"age range: {ages.min_age.min()}-{ages.max_age.max()} BP "
      f"({1950-ages.min_age.min():.0f} to {1950-ages.max_age.max():.0f} BCE)")
