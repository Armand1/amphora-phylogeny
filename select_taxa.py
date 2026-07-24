import numpy as np, pandas as pd

BASE="/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"
z=np.load(f"{BASE}/datafiles/karcher_distmat.npz", allow_pickle=True)
tn=[str(x) for x in z['type_name']]; Da=z['Da']
idx={n:i for i,n in enumerate(tn)}

md=pd.read_csv(f"{BASE}/2026_amphora_297_types_profiles_with_metadata.csv")
md=md.drop_duplicates('type').set_index('type')
# align
meta=md.loc[tn, ['shape_general','shape','place','date_earliest','date_latest']].reset_index()
meta.rename(columns={'index':'type'}, inplace=True)

def fps(names, k):
    """farthest-point sampling within a subset, seeded by most-distant pair."""
    ids=[idx[n] for n in names]
    sub=Da[np.ix_(ids,ids)]
    if len(ids)<=k: return names
    # seed with the two farthest-apart
    i,j=np.unravel_index(np.argmax(sub), sub.shape)
    sel=[i,j]
    while len(sel)<k:
        # min dist from each candidate to selected set
        d=sub[:,sel].min(axis=1)
        d[sel]=-1
        sel.append(int(np.argmax(d)))
    return [names[s] for s in sel]

def pick(mask, k, label):
    names=meta.loc[mask,'type'].tolist()
    chosen=fps(names, k)
    print(f"  {label:28s} pool {len(names):3d} -> {len(chosen)}")
    return chosen

print("== per-class farthest-point subsample (6 greek + 6 italiote) ==")
sel={}
g=meta['place']=="greek"; it=meta['place']=="greek_italiote"

sel['pelike_greek']    = pick((meta.shape_general=="pelike")&g, 6, "pelike greek")
sel['pelike_ital']     = pick((meta.shape_general=="pelike")&it,6, "pelike italiote")
sel['loutro_greek']    = pick((meta.shape_general=="loutrophoros")&g, 6, "loutrophoros greek")
sel['loutro_ital']     = pick((meta.shape_general=="loutrophoros")&it,6, "loutrophoros italiote")
sel['pan_greek']       = pick((meta["shape"]=="amphora_panathenaic")&g, 6, "panathenaic greek")
sel['pan_ital']        = pick((meta["shape"]=="amphora_panathenaic")&it,6, "panathenaic italiote")

# outgroup: 3 deep geometric mainland neck amphorae (oldest), spread by FPS
deep=meta[(meta["shape"]=="amphora_neck")&g&(meta.date_earliest<=-1000)]
sel['outgroup']=fps(deep['type'].tolist(), 3)
print(f"  {'outgroup (deep geometric)':28s} pool {len(deep):3d} -> {len(sel['outgroup'])}")

allsel=[t for v in sel.values() for t in v]
print(f"\nTOTAL TIPS: {len(allsel)}")

out=meta.set_index('type').loc[allsel].reset_index()
out['group']=sum([[k]*len(v) for k,v in sel.items()],[])
out['role']=out['group'].apply(lambda s:'outgroup' if s=='outgroup' else 'test')
out=out[['type','group','role','shape_general','shape','place','date_earliest','date_latest']]
out.to_csv("/private/tmp/claude-501/-Users-armandleroi/b896ff2f-a814-4442-a7e6-244920fe145a/scratchpad/selected_taxa.csv", index=False)
print("\n== selected taxa ==")
pd.set_option('display.max_rows',60,'display.width',160)
print(out.to_string(index=False))
