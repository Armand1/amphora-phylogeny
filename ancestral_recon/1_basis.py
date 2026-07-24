import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from fdasrsf.curve_stats import fdacurve
import fdasrsf.curve_functions as cf
BASE="/Users/armandleroi/Documents/Amphora_phylogeny"
TREE=f"{BASE}/phylogenies/2026_07_three_clade_tree"; DATA=f"{BASE}/Making amphora types/05_2026_amphora_297_types"
sel=pd.read_csv(f"{TREE}/selected_taxa.csv"); prof=pd.read_csv(f"{DATA}/2026_amphora_297_types_profiles_with_metadata.csv")
types=list(sel.type)
beta=np.stack([prof[prof.type==t].sort_values("point_order")[["x","y"]].to_numpy().T for t in types],axis=2)
obj=fdacurve(beta, mode='O', N=beta.shape[1], scale=False)
obj.karcher_mean(rotation=False); obj.srvf_align(rotation=False); obj.shape_pca(no=38)
U,coef,qm,betan = obj.U, obj.coef, obj.q_mean, obj.betan
N=qm.shape[1]
np.savez(f"{TREE}/prep/srvf_recon_basis.npz", U=U, q_mean=qm, coef=coef, types=np.array(types))

def scores_to_curve(c):
    q=qm+(U@c).reshape(2,N); return cf.q_to_curve(q)

# proper overlay: reconstructed vs the ALIGNED normalised curve betan (what scores encode)
fig,ax=plt.subplots(1,4,figsize=(12,4))
for k,i in enumerate([0,10,20,30]):
    tg=betan[:,:,i]; cr=scores_to_curve(coef[:,i])
    tg=tg-tg.mean(1,keepdims=True); cr=cr-cr.mean(1,keepdims=True)
    ax[k].plot(tg[0],tg[1],'k-',lw=2,label='aligned target'); ax[k].plot(cr[0],cr[1],'r--',lw=1.5,label='from scores')
    ax[k].set_title(types[i].split('_')[1][:9]); ax[k].axis('equal'); ax[k].axis('off')
ax[0].legend(fontsize=8); plt.tight_layout(); plt.savefig("srvf_roundtrip_betan.png",dpi=110)
print("saved basis -> prep/srvf_recon_basis.npz ; overlay -> srvf_roundtrip_betan.png")
