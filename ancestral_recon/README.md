# Ancestral profile reconstruction (SRVF)

Draw reconstructed ancestral vase silhouettes at internal tree nodes. Proves SRVF
IS back-projectable: score -> tangent vector (U@c) -> exp map (q_mean + v) ->
inverse-SRVF (q_to_curve) -> curve. Tip round-trip validated (silhouettes overlay
originals; prep/srvf_recon_basis.npz).

Pipeline (SRVF-38 tree, 39 taxa):
1. `1_basis.py`      -> prep/srvf_recon_basis.npz (U, q_mean, coef) + tip round-trip check
   also: prep/srvf38_tip_scores.csv
2. `2_asr.R`         BM ASR (ape::ace, pic) of the 38 scores -> prep/srvf38_node_scores.csv
3. `3_reconstruct.py` node scores + basis -> prep/srvf38_node_curves.csv (height-normalised)
4. `4_draw.R`        ggtree + ancestral silhouettes -> runs/srvf38pc/plots/tree_ancestral_profiles.png

Caveats: reconstruction is a tangent-space linear approx (degrades far from the
Karcher mean); size is discarded (SRVF scale=False) so profiles are height-normalised.
SRVF reconstructions stay on the shape manifold (unlike naive coord-PCA, which can
emit self-intersecting curves at extreme scores).
