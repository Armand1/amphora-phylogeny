import numpy as np, pandas as pd, fdasrsf.curve_functions as cf
B="/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree"
z=np.load(f"{B}/prep/srvf_recon_basis.npz",allow_pickle=True); U,qm=z['U'],z['q_mean']; N=qm.shape[1]
ns=pd.read_csv(f"{B}/prep/srvf38_node_scores.csv",index_col=0)
def recon(c):
    q=qm+(U@c).reshape(2,N); b=cf.q_to_curve(q); b=b.copy()
    b[0]-=b[0].mean(); b/=(b[1].max()-b[1].min()); b[1]-=b[1].min(); return b
rows=[dict(node=int(n),point_order=k,x=b[0,k],y=b[1,k])
      for n,c in ns.iterrows() for b in [recon(c.values.astype(float))] for k in range(b.shape[1])]
pd.DataFrame(rows).to_csv(f"{B}/prep/srvf38_node_curves.csv",index=False)
print("reconstructed", ns.shape[0], "ancestral curves")
