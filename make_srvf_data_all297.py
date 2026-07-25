#!/usr/bin/env python3
"""make_srvf_data_all297.py -- SRVF tangent-PCA for ALL 297 amphora types.
Computes the tangent-PCA once (76 components) and writes BOTH the 38-PC and 76-PC
NEXUS files (38 = first 38 columns of the same basis) + shared taxa TSV.
Fixes the 2 zero-date-range types (spread +/-5 yr) that break RevBayes' tip prior.
Run: python3 make_srvf_data_all297.py
"""
import os, numpy as np, pandas as pd
from fdasrsf.curve_stats import fdacurve
HERE="/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
DATA="/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"
sel=pd.read_csv(f"{HERE}/selected_taxa_all297.csv")
prof=pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv")
types=list(sel.type)
beta=np.stack([prof[prof.type==t].sort_values("point_order")[["x","y"]].to_numpy().T for t in types],axis=2)
print(f"SRVF Karcher mean + tangent-PCA on {len(types)} curves (a few min)...")
obj=fdacurve(beta, mode='O', N=beta.shape[1], scale=False)
obj.karcher_mean(rotation=False); obj.srvf_align(rotation=False); obj.shape_pca(no=76)
S=np.asarray(obj.coef).T[:, :76]                      # (297, 76)
assert S.shape[0]==len(types)
ve=np.asarray(obj.s); ve=ve/ve.sum()
print(f"PC1 {100*ve[0]:.1f}%  | cum@38 {100*np.cumsum(ve)[37]:.1f}%  | cum@76 {100*np.cumsum(ve)[75]:.1f}%")

def write_nex(npc):
    nex=f"{HERE}/data/mini_character_all297_srvf_{npc}.nex"; pad=max(len(t) for t in types)+4
    with open(nex,"w") as f:
        f.write(f"#NEXUS\n[SRVF tPCA all-297, {npc} PCs, make_srvf_data_all297.py]\nBEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={len(types)} NCHAR={npc};\n")
        f.write("  FORMAT DATATYPE=CONTINUOUS MISSING=? GAP=- INTERLEAVE=NO;\n  MATRIX\n")
        for t,row in zip(types,S[:,:npc]):
            f.write("    "+t.ljust(pad)+"\t".join(f"{v:.12f}" for v in row)+"\n")
        f.write("  ;\nEND;\n")
    assert not any("e" in f"{v:.12f}" for v in S[:,:npc].flatten())
    print(f"  wrote {os.path.basename(nex)} (NTAX={len(types)} NCHAR={npc})")
os.makedirs(f"{HERE}/data",exist_ok=True)
write_nex(38); write_nex(76)

ages=prof[prof.type.isin(types)].drop_duplicates("type")[["type","date_earliest","date_latest"]]
ages=ages.assign(min_age=1950-ages.date_latest,max_age=1950-ages.date_earliest)
# fix zero-range types: min_age==max_age degenerates RevBayes' dnUniform tip-age
# prior. Spread by +/-5 yr (range 10). (2 types: 227 @ -175, 272 @ -250.)
zr=ages.max_age==ages.min_age
ages.loc[zr,"min_age"]-=5; ages.loc[zr,"max_age"]+=5
print(f"  zero-range types spread to 10yr: {int(zr.sum())}")
ages[["type","min_age","max_age"]].rename(columns={"type":"taxon"}).to_csv(
    f"{HERE}/data/mini_taxa_all297.tsv",sep="\t",index=False)
print(f"  wrote mini_taxa_all297.tsv ({len(ages)} taxa); age {ages.min_age.min():.0f}-{ages.max_age.max():.0f} BP")
