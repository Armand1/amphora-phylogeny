# Three-clade amphora test tree

A purpose-built RevBayes tree to test the **Athens → Italy rate hypothesis**:
when Athenian potters moved to Italy, did every vessel class undergo rapid shape
change? Three distinctive classes, each present on both sides of the divide, give
three independent replicates of "tradition transplanted."

## Taxa (39 tips)

| clade | Greek | Italiote | role |
|---|---|---|---|
| pelike | 6 | 6 | test |
| loutrophoros | 6 | 6 | test |
| panathenaic (amphora) | 6 | 6 | test |
| deep geometric neck amphorae | 3 | — | outgroup (roots tree + polarises geography) |

Selected by farthest-point sampling within each place×class subgroup on the
Karcher-SRVF distance matrix (`select_taxa.py` → `selected_taxa.csv`). Ages span
−1040 to −100 BCE; classes overlap in time across the divide, so geography is not
confounded with age.

## Quick NJ look (before the dated model)

`three_clade_nj.Rmd` → `three_clade_nj.html`. PCA of the 39 aligned outlines:
PC1 = 67.6% (squatness), 10 PCs = 98.9%, 20 PCs = 99.9% (saturated). On a naive
NJ tree **no class is monophyletic** at any PC count (5/10/20), and the topology
is fully stable by 20 PCs. Loutrophoros comes closest (4 interlopers);
panathenaics and pelikes are dispersed. Conclusion: model-free shape distance
cannot recover the classes — the test is whether **FBD + tip-dates + BM** can.

## RevBayes files

- `data/mini_character.nex` — 39 tips × 10 continuous characters (coordinate-PCA
  scores of the aligned outlines). Regenerate with a different PC count via
  `Rscript make_revbayes_data.R 20`.
- `data/mini_taxa.tsv` — min/max tip ages in years BP, datum 1950
  (`min_age = 1950 − date_latest`, `max_age = 1950 − date_earliest`).
- `master_script.Rev` — constant-rate FBD + strict clock + BM on the 10 PCs.
  Edit `BASE` for the run location. Produces `output/simple_run*.{log,trees}`
  and `output/simple_run_mcc.tree`.
- `scripts/FBD.Rev`, `strict_clock.Rev`, `simple_BM.Rev` — model modules.

### Changes from the pilot's model
- `trait <- 1:10` (this tree has 10 PCs, not 65).
- **`origin_time` ceiling raised 3100 → 4000 BP.** In the 65-PC pilot the root
  jammed against the 3100 cap and the chronology broke; our oldest tip is 2990 BP.
- `psi` / `lambda` priors left at the pilot values but flagged in `FBD.Rev`
  (suspected Myr-scale priors on year-scale data — raise with Joel).

### Character basis caveat
The characters are **coordinate-PCA** of the Karcher-aligned outlines, computed on
just these 39 types — self-consistent with the NJ analysis but a *different basis*
from the pilot's SRVF-tPCA. Fine for this standalone test; if we later want strict
cross-registration with the pilot or other trees, recompute through the SRVF-tPCA
pipeline.

## Workflow (as for the pilot)
1. Run `master_script.Rev` in RevBayes → constant-rate MCC tree + posterior trees.
2. Diagnostics in Tracer (ESS/PSRF, watch `origin_time`, `psi`, `sigma2`).
3. Apply relaxed (UCLN) branch rates on the fixed topology (separate step).
4. Tempo pipeline: geography ASR → classify intervals → transplant-branch rate
   test across the three clades (see `../Tempo methodology …md`).
