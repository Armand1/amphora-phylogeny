import numpy as np, pandas as pd
BASE="/Users/armandleroi/Documents/Amphora_phylogeny"
TREE=f"{BASE}/phylogenies/2026_07_three_clade_tree"
DATA=f"{BASE}/Making amphora types/05_2026_amphora_297_types"
z=np.load(f"{DATA}/datafiles/karcher_distmat.npz",allow_pickle=True)
tn=[str(x) for x in z['type_name']]; Da=z['Da']; idx={n:i for i,n in enumerate(tn)}
md=pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv").drop_duplicates('type').set_index('type')

def fps(names,k):
    ids=[idx[n] for n in names]; sub=Da[np.ix_(ids,ids)]
    if len(ids)<=k: return names
    i,j=np.unravel_index(np.argmax(sub),sub.shape); sel=[i,j]
    while len(sel)<k:
        d=sub[:,sel].min(axis=1); d[sel]=-1; sel.append(int(np.argmax(d)))
    return [names[s] for s in sel]

rows=[]
def add(t,group,role):
    r=md.loc[t]; rows.append(dict(type=t,group=group,role=role,shape_general=r.shape_general,
        shape=r['shape'],place=r.place,date_earliest=r.date_earliest,date_latest=r.date_latest))

# ALL neck + ALL belly (zero selection)
for t in md.index:
    if md.loc[t,'shape']=='amphora_neck':  add(t,'neck_greek' if md.loc[t,'place']=='greek' else 'neck_ital','test')
    if md.loc[t,'shape']=='amphora_belly': add(t,'belly_greek' if md.loc[t,'place']=='greek' else 'belly_ital','test')

# 12 anchors: 4 each pelike/pan/loutro by farthest-point (fixed reference set)
def anchors(mask,tag,k=4):
    names=[t for t in md.index if mask(t)]
    for t in fps(names,k): add(t,f'anchor_{tag}','anchor')
anchors(lambda t: md.loc[t,'shape_general']=='pelike','pelike')
anchors(lambda t: md.loc[t,'shape']=='amphora_panathenaic','panathenaic')
anchors(lambda t: md.loc[t,'shape_general']=='loutrophoros','loutrophoros')

nb=pd.DataFrame(rows).drop_duplicates('type')
nb.to_csv(f"{TREE}/selected_taxa_neckbelly.csv",index=False)
nb['place2']=nb.place.map({'greek':'Greek','greek_italiote':'Italiote'})
nb['cls']=nb.group.str.replace('_greek','').str.replace('_ital','').str.replace('anchor_','')
print("=== neck+belly tree composition ===")
print(pd.crosstab(nb.cls, nb.place2, margins=True))
print(f"\nTOTAL: {len(nb)} taxa  (test={sum(nb.role=='test')}, anchors={sum(nb.role=='anchor')})")
