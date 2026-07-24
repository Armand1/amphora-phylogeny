# Three-clade amphora test tree

RevBayes FBD tree testing the **Athens → Italy rate hypothesis**: three
distinctive vessel classes (pelike / loutrophoros / panathenaic), each present on
both sides of the divide, rooted on deep geometric neck amphorae. 39 tips.

This folder is a git repo. On the HPC it is cloned as `amphora-phylogeny` and the
PBS scripts reference paths like `amphora-phylogeny/master_script.Rev` — so the
**model-source files must stay at the repo root** (do not move them into subfolders).

## Layout

```
├── master_script.Rev        top-level model (edit `trait <- 1:N` for PC count)
├── mcc_only.Rev             rebuild MCC tree if the summary step was interrupted
├── submit_hpc.pbs           HPC job: the full MCMC
├── submit_mcc.pbs           HPC job: MCC-tree rebuild only (single process)
├── scripts/                 FBD.Rev · strict_clock.Rev · simple_BM.Rev
├── data/                    mini_character.nex (39 × N PCs) · mini_taxa.tsv (ages, BP)
├── output/                  HPC writes here at runtime (git-empty)
│
├── make_revbayes_data.R     profiles → data/*.nex + *.tsv  (arg = #PCs, default 10)
├── select_taxa.py           farthest-point taxon selection → selected_taxa.csv
├── selected_taxa.csv        the 39 chosen types
│
├── three_clade_nj.Rmd       pre-analysis: PCA + naive NJ tree (PC-independent)
├── diagnostics.Rmd          Tracer-style MCMC diagnostics for a run
├── tree_plots.Rmd           dated tree with pot profiles + HPD bars
│
├── prep/                    one-off prep artifacts (NJ html, sub distance matrix)
└── runs/                    RESULTS, one folder per run (local only, git-ignored)
    ├── 10pc/  output/ · plots/ · *.html
    └── 20pc/  output/ · plots/
```

`runs/` and all `*.html` are git-ignored: they are regenerated, not source. The
git repo carries only the model + the notebooks that produce the results.

## Workflow

1. **Build data** (choose PC count; also set `trait <- 1:N` in `master_script.Rev`):
   ```bash
   Rscript make_revbayes_data.R 20
   git add -A && git commit -m "20 PCs" && git push
   ```
2. **Run on HPC** (see `../../notes/HPC_idiots_guide.md` for the full walk-through):
   ```bash
   cd ~/GitHub/amphora-phylogeny && git pull
   cd ~/GitHub && qsub amphora-phylogeny/submit_hpc.pbs
   qsub amphora-phylogeny/submit_mcc.pbs      # after it finishes, rebuilds the MCC tree
   ```
3. **Pull results back** into the matching run folder (run in a Mac terminal):
   ```bash
   scp -r amleroi@login.cx3.hpc.imperial.ac.uk:'~/GitHub/amphora-phylogeny/output' \
         runs/20pc/output
   ```
4. **Analyse**: render `diagnostics.Rmd` and `tree_plots.Rmd` with `run_dir` set to
   `runs/<pc>/output` (output the html/plots into that run folder).

## Status (2026-07-24)

- **10-PC run complete.** Deepest split (posterior 1.00) = squat (pelike) vs
  elongate (loutrophoros + panathenaic merged); Greek & Italiote intermingle
  within every clade. Diagnostics clean; sigma2 off the floor; root ≈ 1281 BCE.
  Within-class nodes < 0.5 support (near-identical shapes → genuinely unresolvable).
- **20-PC run** in progress, to test whether more characters buy within-class
  resolution or whether it is a hard limit of the data.
