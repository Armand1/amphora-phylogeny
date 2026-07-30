#!/bin/sh
# Pull the good238 tuned run results down from cx3 into runs/good238_{38,76}_tuned/output/.
#
# Deliberately does NOT fetch the six *_run_N.trees posteriors (~175 MB each,
# ~1 GB total). They stay on the cluster: the MCC tree is built HPC-side by
# submit_mcc_good238_*_tuned_perrun.pbs, and diagnostics_tuned.Rmd needs only
# the .log files. Run with --trees if you want the posteriors locally as well
# (needed for multi-tree tempo analysis, not for diagnostics or tree plots).
#
# Usage:   sh pull_good238_tuned.sh [--trees]
# Run this ON THE MAC, with the VPN up.

set -e
HPC="amleroi@login.cx3.hpc.imperial.ac.uk"
REMOTE="~/GitHub/amphora-phylogeny/output"
HERE="$(cd "$(dirname "$0")" && pwd)"

WANT_TREES=0
[ "$1" = "--trees" ] && WANT_TREES=1

for pc in 38 76; do
  RUN="good238_${pc}_tuned"
  DEST="$HERE/runs/$RUN/output"
  mkdir -p "$DEST"
  echo "=== $RUN -> runs/$RUN/output"

  # per-run parameter logs (small) — what diagnostics_tuned.Rmd reads
  scp "$HPC:$REMOTE/simple_run_${RUN}_run_*.log" "$DEST/" || echo "  (no run logs)"

  # salvage MCC tree + its build log, and the clean-finish MCC if it exists
  scp "$HPC:$REMOTE/simple_run_${RUN}_mcc_perrun.tree" "$DEST/" || echo "  (no perrun MCC yet)"
  scp "$HPC:$REMOTE/mcc_logfile_${RUN}_perrun.log"     "$DEST/" || echo "  (no perrun MCC log)"
  scp "$HPC:$REMOTE/simple_run_${RUN}_mcc.tree"        "$DEST/" 2>/dev/null || true

  # the MCMC console log (errors, walltime kill message)
  scp "$HPC:$REMOTE/logfile_${RUN}.log" "$DEST/" || echo "  (no logfile)"

  if [ "$WANT_TREES" -eq 1 ]; then
    echo "  fetching posteriors (~525 MB for this run)..."
    scp "$HPC:$REMOTE/simple_run_${RUN}_run_*.trees" "$DEST/"
  fi
done

echo
echo "Done. Local sizes:"
du -sh "$HERE"/runs/good238_*_tuned/output 2>/dev/null || true
echo
echo "Next: render diagnostics_tuned.Rmd for each run, e.g."
echo "  Rscript -e 'rmarkdown::render(\"diagnostics_tuned.Rmd\","
echo "    params=list(run_dir=\"runs/good238_38_tuned/output\","
echo "                run_prefix=\"simple_run_good238_38_tuned\"),"
echo "    output_file=\"runs/good238_38_tuned/diagnostics_38_tuned.html\")'"
