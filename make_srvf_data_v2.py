#!/usr/bin/env python3
"""make_srvf_data_v2.py -- SRVF tangent-PCA for the v2 Greek/Italiote dataset.

330 taxa = 323 types (amphora 199 + hydria 95 + stamnos 29) + 7 SOS vases,
built at 'Making amphora types/08_2026_greek_italiote_v2/'.  Sibling of
make_srvf_data_good238.py; that script and its data files are untouched.

The tangent-PCA basis is recomputed on THIS taxon set, so its PC axes are not
the good-238 axes -- the two trees are not comparable coordinate-for-coordinate.

The 7 SOS vases are individual pots, not Karcher-mean types, and enter as a
fixed taxon set (they must never be run through density reduction).  All curves
are 250-point open one-sided outlines from the same pipeline, so they stack
directly.

Prints the cumulative-variance curve so the PC count can be chosen rather than
inherited: 38 was right for good-238, but this dataset has four shape classes
where that had three, and variance spread across more genuine shape modes means
38 PCs need not reach the same coverage.

Run: python3 make_srvf_data_v2.py [npc ...]      (default: report only)
"""
import os
import sys

import numpy as np
import pandas as pd
from fdasrsf.curve_stats import fdacurve

HERE = "/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
MT = "/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types"
V2 = f"{MT}/08_2026_greek_italiote_v2"
PROF = "/Users/armandleroi/Documents/Global-Pot-Project/main_data/all_profiles_310726.parquet"

KARCHER = {
    'amphora': f"{MT}/05_2026_amphora_297_types/datafiles/types_karcher_profiles.csv",
    'hydria':  f"{MT}/07_2026_hydria_types/types_karcher_profiles_hydria.csv",
    'stamnos': f"{MT}/07_2026_stamnos_types/types_karcher_profiles_stamnos.csv",
}

# ---------------------------------------------------------------- assemble
keep = pd.read_csv(f"{V2}/type_list_v2.csv", keep_default_na=False)
sos = pd.read_csv(f"{V2}/sos_taxa.csv", keep_default_na=False)

names, curves, cls, d_early, d_late = [], [], [], [], []

for cl, path in KARCHER.items():
    k = pd.read_csv(path)
    sub = keep[keep.type_class == cl]
    have = set(k.type_global_id)
    missing = set(sub.type_global_id) - have
    assert not missing, f"{cl}: no Karcher profile for {sorted(missing)}"
    for r in sub.sort_values('type_global_id').itertuples():
        c = k[k.type_global_id == r.type_global_id].sort_values('point_order')
        assert len(c) == 250, f"{r.type} has {len(c)} points"
        names.append(r.type)
        curves.append(c[['x', 'y']].to_numpy().T)
        cls.append(cl)
        d_early.append(r.date_earliest)
        d_late.append(r.date_latest)

p = pd.read_parquet(PROF)
p = p[p.gpp_no.isin(sos.gpp_no)]
for r in sos.sort_values('mid', ascending=False).itertuples():
    c = p[p.gpp_no == r.gpp_no].sort_values('point_order')
    assert len(c) == 250, f"{r.gpp_no} has {len(c)} points"
    # tip label follows the type convention so downstream parsers still work
    names.append(f"greek_sos_iron_age_{int(r.mid)}BCE_{r.gpp_no}")
    curves.append(c[['x', 'y']].to_numpy().T)
    cls.append('sos')
    d_early.append(float(r.date_earliest))
    d_late.append(float(r.date_latest))

assert len(names) == len(set(names)), "duplicate taxon names"
beta = np.stack(curves, axis=2)
print(f"assembled {beta.shape[2]} taxa  "
      f"{pd.Series(cls).value_counts().to_dict()}")

# ---------------------------------------------------------------- tangent PCA
NMAX = 100
print(f"SRVF Karcher mean + tangent-PCA on {len(names)} curves (several min)...")
obj = fdacurve(beta, mode='O', N=beta.shape[1], scale=False)
obj.karcher_mean(rotation=False)
obj.srvf_align(rotation=False)
obj.shape_pca(no=NMAX)
S = np.asarray(obj.coef).T[:, :NMAX]
assert S.shape[0] == len(names), f"coef shape {S.shape} vs {len(names)} taxa"

ve = np.asarray(obj.s)
ve = ve / ve.sum()
cum = np.cumsum(ve)
print(f"\nPC1 {100*ve[0]:.1f}%")
print("cumulative variance:")
for n in (5, 10, 20, 30, 38, 50, 60, 76, 100):
    if n <= len(cum):
        print(f"  {n:3d} PCs  {100*cum[n-1]:6.2f}%   (PC{n} alone {100*ve[n-1]:.3f}%)")
for thr in (0.95, 0.975, 0.978, 0.99):
    print(f"  {100*thr:.1f}% reached at {int(np.searchsorted(cum, thr)) + 1} PCs")
np.savez(f"{HERE}/prep/srvf_v2_scores.npz", scores=S, taxon=np.array(names),
         type_class=np.array(cls), var_explained=ve)
print(f"\nwrote prep/srvf_v2_scores.npz")

# ---------------------------------------------------------------- outputs
def write_nex(npc):
    nex = f"{HERE}/data/mini_character_v2_srvf_{npc}.nex"
    pad = max(len(t) for t in names) + 4
    with open(nex, "w") as f:
        f.write(f"#NEXUS\n[SRVF tPCA v2 Greek/Italiote, {npc} PCs, make_srvf_data_v2.py]\n")
        f.write("BEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={len(names)} NCHAR={npc};\n")
        f.write("  FORMAT DATATYPE=CONTINUOUS MISSING=? GAP=- INTERLEAVE=NO;\n  MATRIX\n")
        for t, row in zip(names, S[:, :npc]):
            f.write("    " + t.ljust(pad) + "\t".join(f"{v:.12f}" for v in row) + "\n")
        f.write("  ;\nEND;\n")
    # RevBayes' readContinuousCharacterData cannot parse scientific notation
    assert not any("e" in f"{v:.12f}" for v in S[:, :npc].flatten())
    print(f"  wrote {os.path.basename(nex)} (NTAX={len(names)} NCHAR={npc})")


ages = pd.DataFrame({'taxon': names,
                     'min_age': [1950 - x for x in d_late],
                     'max_age': [1950 - x for x in d_early]})
zr = ages.max_age == ages.min_age          # zero range breaks the tip-age prior
ages.loc[zr, 'min_age'] -= 5
ages.loc[zr, 'max_age'] += 5
print(f"\nzero-range taxa spread to 10yr: {int(zr.sum())}")

if len(sys.argv) > 1:
    os.makedirs(f"{HERE}/data", exist_ok=True)
    for a in sys.argv[1:]:
        write_nex(int(a))
    ages.to_csv(f"{HERE}/data/mini_taxa_v2.tsv", sep="\t", index=False)
    print(f"  wrote mini_taxa_v2.tsv ({len(ages)} taxa); "
          f"age {ages.min_age.min():.0f}-{ages.max_age.max():.0f} BP")
else:
    print("(report only -- pass PC counts as arguments to write NEXUS + taxa TSV)")
