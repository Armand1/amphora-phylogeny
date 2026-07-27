#!/usr/bin/env python3
"""make_srvf_data_good238.py -- SRVF tangent-PCA for the 238 GOOD amphora types.
Post-QC dataset: 46 'bad' types (chopped lips/bases) and 13 'ok but not useful'
(coarseware/domestic) types removed by eye. Recomputes the tangent-PCA basis on
the survivors only (so the PC axes no longer carry the junk), writes 38-PC and
76-PC NEXUS + taxa TSV. Sibling of make_srvf_data_all297.py.
Run: python3 make_srvf_data_good238.py
"""
import os, numpy as np, pandas as pd
from fdasrsf.curve_stats import fdacurve
HERE="/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
DATA="/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"
prof=pd.read_csv(f"{DATA}/2026_amphora_good_types_profiles_with_metadata.csv")
types=prof.drop_duplicates("type").sort_values("type_global_id").type.tolist()
beta=np.stack([prof[prof.type==t].sort_values("point_order")[["x","y"]].to_numpy().T for t in types],axis=2)
print(f"SRVF Karcher mean + tangent-PCA on {len(types)} good curves (a few min)...")
obj=fdacurve(beta, mode='O', N=beta.shape[1], scale=False)
obj.karcher_mean(rotation=False); obj.srvf_align(rotation=False); obj.shape_pca(no=76)
S=np.asarray(obj.coef).T[:, :76]
assert S.shape[0]==len(types)
ve=np.asarray(obj.s); ve=ve/ve.sum()
print(f"PC1 {100*ve[0]:.1f}%  | cum@38 {100*np.cumsum(ve)[37]:.1f}%  | cum@76 {100*np.cumsum(ve)[75]:.1f}%")

def write_nex(npc):
    nex=f"{HERE}/data/mini_character_good238_srvf_{npc}.nex"; pad=max(len(t) for t in types)+4
    with open(nex,"w") as f:
        f.write(f"#NEXUS\n[SRVF tPCA good-238, {npc} PCs, make_srvf_data_good238.py]\nBEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={len(types)} NCHAR={npc};\n")
        f.write("  FORMAT DATATYPE=CONTINUOUS MISSING=? GAP=- INTERLEAVE=NO;\n  MATRIX\n")
        for t,row in zip(types,S[:,:npc]):
            f.write("    "+t.ljust(pad)+"\t".join(f"{v:.12f}" for v in row)+"\n")
        f.write("  ;\nEND;\n")
    assert not any("e" in f"{v:.12f}" for v in S[:,:npc].flatten())
    print(f"  wrote {os.path.basename(nex)} (NTAX={len(types)} NCHAR={npc})")
os.makedirs(f"{HERE}/data",exist_ok=True)
write_nex(38); write_nex(76)

ages=prof.drop_duplicates("type")[["type","date_earliest","date_latest"]]
ages=ages.assign(min_age=1950-ages.date_latest,max_age=1950-ages.date_earliest)
zr=ages.max_age==ages.min_age                     # zero-range breaks the tip-age prior
ages.loc[zr,"min_age"]-=5; ages.loc[zr,"max_age"]+=5
print(f"  zero-range types spread to 10yr: {int(zr.sum())}")
ages[["type","min_age","max_age"]].rename(columns={"type":"taxon"}).to_csv(
    f"{HERE}/data/mini_taxa_good238.tsv",sep="\t",index=False)
print(f"  wrote mini_taxa_good238.tsv ({len(ages)} taxa); age {ages.min_age.min():.0f}-{ages.max_age.max():.0f} BP")
