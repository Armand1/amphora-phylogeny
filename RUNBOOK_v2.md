# v2 Greek/Italiote tree — runbook

330 taxa. Dataset built at
`Making amphora types/08_2026_greek_italiote_v2/` (see its README for the
exclusion rules). All files here are **siblings** — nothing belonging to the
good-238 or earlier configs was modified.

## What the 330 tips actually are

`type_class` in the dataset names the *build*, not the vessel shape. The
"amphora" build has always carried three nominal shapes. Across the tips:

| nominal shape | n |
|---|---|
| amphora | 107 |
| hydria | 95 |
| pelike | 67 |
| stamnos | 29 |
| loutrophoros | 25 |
| SOS | 7 |

## Characters

`make_srvf_data_v2.py` recomputes the SRVF tangent-PCA basis on this taxon set.
**The PC axes are therefore not the good-238 axes** — the two trees cannot be
compared coordinate-for-coordinate.

```
PC1 27.6%   |   38 PCs 97.11%   |   45 PCs 97.80%   |   76 PCs 99.45%
```

Two arms, run concurrently:

- **38 PC** — same character count as the good-238 reference tree, so a
  difference is attributable to the dataset rather than to dimensionality.
- **45 PC** — matched *variance coverage* to good-238 (97.8%). 38 was never a
  magic number; it was the largest count that kept sigma² off its old floor.
  With the floor lowered that constraint is gone, and this arm tests whether
  the extra coverage is now affordable.

Note PC1 rose (23.7% → 27.6%) even though the dataset gained shape classes. The
extra classes lengthened the tail rather than flattening the head.

## Model changes vs good-238

- `scripts/strict_clock_v2.Rev` — `sigma2 ~ dnLoguniform(1e-6, 1)`, was 1e-3.
  The good-238 38-PC run cleared the old floor by only 12% (min 0.00112), and
  on a new basis sigma²'s scale shifts. Not prior-neutral — see the comment in
  the file.
- `FBD_tuned.Rev` and `simple_BM.Rev` unchanged.
- `printgen` 200 → 400, so 4M generations yields the same ~10k samples/run that
  good-238 produced at 2M — a tree-file size already known to summarise.
- mem 32 → 64 GB.

## Checkpointing — the point of it

`mymcmcmc.run(generations=4000000, checkpointInterval=20000, checkpointFile=…)`.

The generation target is deliberately beyond any one wall. **ASDSF decides when
to stop, not the wall.** 2M generations was calibrated for 238 tips; 330 tips
explore a larger tree space and need more, not the same. Run to the wall,
resubmit the `_restart` config, repeat.

On restart the monitors switch to `append=TRUE` — without that a resumed job
silently truncates everything the previous job wrote.

**Nothing in this repo has ever used RevBayes checkpointing, and mcmcmc
checkpoint/restart is fussier than plain mcmc** (heated-chain state, per-run
files). Prove it on the smoke test before committing days of compute.

## Order of operations

```
# 0. push, then on the HPC:  cd ~/GitHub/amphora-phylogeny && git pull

# 1. SMOKE TEST checkpointing (30 min wall, 3000 gens on the real data)
qsub submit_hpc_v2_ckptest.pbs
qsub submit_hpc_v2_ckptest_restart.pbs
#    PASS if the restart log continues past generation 3000 rather than
#    starting at 1, and simple_run_v2_ckptest.log GREW rather than being
#    truncated. If it fails, fix it here — not at 72 h.

# 2. the real runs, concurrently
qsub submit_hpc_v2_38.pbs
qsub submit_hpc_v2_45.pbs

# 3. after each walltime kill, resubmit as often as needed
qsub submit_hpc_v2_38_restart.pbs
qsub submit_hpc_v2_45_restart.pbs

# 4. summarise (per-run: a checkpointed run stopped at the wall by design
#    never writes the combined .trees)
qsub submit_mcc_v2_38_perrun.pbs
qsub submit_mcc_v2_45_perrun.pbs

# 5. pull, then ASDSF locally
python3 compute_asdsf.py --thin 10 …
```

## Sizing expectation

good-238 managed ~26,000 gens/h at 238 tips. The BM pruning likelihood is
O(tips × characters) and topology moves cost more on a bigger tree, so expect
≥1.4× slower — roughly 18–19k gens/h, i.e. **~105 h for 2M generations**, more
for the 45-PC arm. A single 72 h wall will not reach 2M. That is expected and
is why the restart configs exist.

If a longer queue than 72 h is available on cx3 it is worth using, but it does
not remove the need for checkpointing — it only reduces the number of
resubmissions.

## Files

| | |
|---|---|
| `make_srvf_data_v2.py` | basis + NEXUS + taxa TSV; run with no args to report the variance curve only |
| `data/mini_character_v2_srvf_{38,45}.nex` | 330 × 38 and 330 × 45 |
| `data/mini_taxa_v2.tsv` | 330 taxa, 2025–3000 BP; 1 zero-range taxon spread to 10 yr |
| `prep/srvf_v2_scores.npz` | scores, taxon names, `type_class`, variance explained |
| `scripts/v2_body.Rev` | shared body; callers set `NPC`, `TAG`, `RESTART` |
| `scripts/strict_clock_v2.Rev` | sigma² floor 1e-6 |
| `master_script_v2_{38,45}[_restart].Rev` | four thin callers |
| `master_script_v2_ckptest[_restart].Rev`, `scripts/v2_body_ckptest.Rev` | smoke test |
| `submit_hpc_v2_*.pbs`, `mcc_only_v2_*_perrun.Rev`, `submit_mcc_v2_*_perrun.pbs` | |
