import numpy as np, pandas as pd
BASE="/Users/armandleroi/Documents/Amphora_phylogeny"
TREE=f"{BASE}/phylogenies/2026_07_three_clade_tree"
DATA=f"{BASE}/Making amphora types/05_2026_amphora_297_types"

z=np.load(f"{DATA}/datafiles/karcher_distmat.npz", allow_pickle=True)
tn=[str(x) for x in z['type_name']]; Da=z['Da']; idx={n:i for i,n in enumerate(tn)}
md=pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv").drop_duplicates('type').set_index('type')
old=pd.read_csv(f"{TREE}/selected_taxa.csv")                 # the current 39
used=set(old.type)

def fps(names,k,seed_used):
    """farthest-point among `names`, seeded to be far from already-chosen seed_used."""
    ids=[idx[n] for n in names]; sub=Da[np.ix_(ids,ids)]
    seed=[i for i,n in enumerate(names) if n in seed_used]
    sel=list(seed) if seed else [int(np.unravel_index(np.argmax(sub),sub.shape)[0])]
    while len(sel)<k+len(seed):
        d=sub[:,sel].min(axis=1); d[sel]=-1; sel.append(int(np.argmax(d)))
    return [names[s] for s in sel if names[s] not in seed_used]   # only the NEW ones

def is_pan(t): return md.loc[t,'shape']=='amphora_panathenaic'
def cls(t):
    if md.loc[t,'shape_general'] in ('pelike','loutrophoros'): return md.loc[t,'shape_general']
    if is_pan(t): return 'panathenaic'
    return md.loc[t,'shape_general']

rows=[]
def add(t, group, role):
    r=md.loc[t]; rows.append(dict(type=t, group=group, role=role, shape_general=r.shape_general,
        shape=r['shape'], place=r.place, date_earliest=r.date_earliest, date_latest=r.date_latest))

# ALL loutrophoroi + ALL panathenaics
for t in md.index:
    if md.loc[t,'shape_general']=='loutrophoros':
        add(t,'loutro_greek' if md.loc[t,'place']=='greek' else 'loutro_ital','test')
    elif is_pan(t):
        add(t,'pan_greek' if md.loc[t,'place']=='greek' else 'pan_ital','test')

# PELIKES: keep the current 12, add 15 (7 Greek + 8 Italiote) by farthest-point
pel=[t for t in md.index if md.loc[t,'shape_general']=='pelike']
pel_g=[t for t in pel if md.loc[t,'place']=='greek']
pel_i=[t for t in pel if md.loc[t,'place']=='greek_italiote']
keep_pel=[t for t in used if t in pel]
newg=fps(pel_g, 6+7, set(keep_pel))[:7]     # 6 kept greek + 7 new
newi=fps(pel_i, 6+8, set(keep_pel))[:8]     # 6 kept ital  + 8 new
for t in keep_pel+newg: add(t,'pelike_greek' if md.loc[t,'place']=='greek' else 'pelike_ital','test')
for t in newi: add(t,'pelike_ital','test')

# OUTGROUP: keep the same 3 deep geometric necks
for t in old[old.role=='outgroup'].type: add(t,'outgroup','outgroup')

big=pd.DataFrame(rows).drop_duplicates('type')
big.to_csv(f"{TREE}/selected_taxa_big.csv", index=False)

big['cls']=big.type.map(cls)
print("=== BIG dataset composition ===")
print(pd.crosstab(big.cls, big.place, margins=True))
print(f"\nTOTAL: {len(big)} taxa")
print(f"supersets the current 39? {used.issubset(set(big.type))}")
