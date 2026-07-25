import pandas as pd
DATA="/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types"
md=pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv").drop_duplicates('type')
def cls(r):
    if r.shape_general in ('pelike','loutrophoros'): return r.shape_general
    if r['shape']=='amphora_panathenaic': return 'panathenaic'
    if r['shape']=='amphora_neck': return 'neck'
    if r['shape']=='amphora_belly': return 'belly'
    return 'amphora_other'
md=md.assign(cls=md.apply(cls,axis=1),
             og=md.place.map({'greek':'greek','greek_italiote':'ital'}))
out=md.assign(group=md.cls+'_'+md.og, role='test')[
    ['type','group','role','shape_general','shape','place','date_earliest','date_latest']]
out.to_csv("selected_taxa_all297.csv",index=False)
out['place2']=out.place.map({'greek':'G','greek_italiote':'I'})
out['cls']=out.group.str.replace('_greek','').str.replace('_ital','')
print("=== ALL 297 taxa by class x place ===")
print(pd.crosstab(out.cls,out.place2,margins=True))
print(f"date NAs: {out.date_earliest.isna().sum()+out.date_latest.isna().sum()}")
