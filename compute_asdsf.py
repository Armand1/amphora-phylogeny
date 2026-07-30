#!/usr/bin/env python3
"""compute_asdsf.py -- topological convergence (ASDSF) across RevBayes runs.

The tuned good238 runs were killed by the walltime, so RevBayes never printed its
end-of-run split-frequency diagnostic. ESS/PSRF on the continuous parameters says
nothing about whether the independent runs agree on TOPOLOGY -- this does.

ASDSF = average, over splits, of the standard deviation of that split's frequency
across the independent runs. Rule of thumb (MrBayes): < 0.01 excellent,
< 0.05 acceptable, > 0.10 not converged.

Dependency-free (stdlib only) so it runs on the HPC without a module environment.
Run it on the cluster, where the .trees files already live -- they are ~175 MB
each and do not need to come down to the Mac.

Usage:
  python3 compute_asdsf.py output/simple_run_good238_38_tuned_run_{1,2,3}.trees
  python3 compute_asdsf.py --burnin 0.25 --thin 5 <files...>
  python3 compute_asdsf.py --selftest
"""
import argparse, math, re, sys
from collections import Counter

ANNOT = re.compile(r"\[[^\]]*\]")          # [&index=..,posterior=..] blocks
NUM   = re.compile(r":\s*-?[0-9.eE+-]+")   # branch lengths


def newick_splits(nwk, taxa_index):
    """Return the set of non-trivial splits in one newick string.

    Each split is a frozenset of taxon indices, canonicalised so it never
    contains taxon 0. Degree-2 nodes (sampled ancestors) collapse naturally
    because we store splits in a set.
    """
    stack = [set()]
    token = []

    def flush():
        if token:
            name = "".join(token).strip().strip("'\"")
            if name:
                stack[-1].add(taxa_index[name])
            token.clear()

    for ch in nwk:
        if ch == "(":
            stack.append(set())
        elif ch == ",":
            flush()
        elif ch == ")":
            flush()
            top = stack.pop()
            stack[-1] |= top
            stack[-1].add(("SPLIT", frozenset(top)))  # marker, stripped below
        else:
            token.append(ch)
    flush()

    # collect the markers we buried in the sets
    out = set()
    n = len(taxa_index)

    def harvest(s):
        for item in s:
            if isinstance(item, tuple):
                out.add(item[1])

    # markers accumulate upward, so the root set holds them all
    harvest(stack[-1])
    # also strip markers from the leaf-index sets we unioned upward
    clean = set()
    for sp in out:
        idx = frozenset(i for i in sp if not isinstance(i, tuple))
        if 1 < len(idx) < n - 1:                      # drop trivial splits
            clean.add(idx if 0 not in idx else frozenset(range(n)) - idx)
    return clean


def extract_newick(line):
    """Pull the newick out of a RevBayes .trees line (tab columns, or NEXUS)."""
    i = line.find("(")
    if i < 0:
        return None
    j = line.rfind(")")
    if j <= i:
        return None
    s = line[i:j + 1]
    s = ANNOT.sub("", s)
    s = NUM.sub("", s)
    return s


def read_taxa(path):
    """Taxon name -> index, from the first tree in the file."""
    with open(path) as fh:
        for line in fh:
            nwk = extract_newick(line)
            if nwk is None:
                continue
            names = [t.strip().strip("'\"") for t in re.split(r"[(),]", nwk)]
            names = [t for t in names if t]
            return {t: i for i, t in enumerate(sorted(set(names)))}
    raise SystemExit(f"no trees found in {path}")


def run_freqs(path, taxa_index, burnin, thin):
    """Split -> frequency, for one run."""
    trees = []
    with open(path) as fh:
        for line in fh:
            nwk = extract_newick(line)
            if nwk is not None:
                trees.append(nwk)
    b = int(len(trees) * burnin)
    trees = trees[b:][::thin]
    if not trees:
        raise SystemExit(f"no post-burnin trees in {path}")
    c = Counter()
    for nwk in trees:
        for sp in newick_splits(nwk, taxa_index):
            c[sp] += 1
    n = len(trees)
    return {sp: k / n for sp, k in c.items()}, n


def asdsf(freqs, minfreq):
    """Average SD of split frequencies across runs, over splits above minfreq."""
    keys = set()
    for f in freqs:
        keys |= {sp for sp, v in f.items() if v >= minfreq}
    if not keys:
        return float("nan"), float("nan"), 0
    sds = []
    for sp in keys:
        vals = [f.get(sp, 0.0) for f in freqs]
        m = sum(vals) / len(vals)
        sds.append(math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals)))
    return sum(sds) / len(sds), max(sds), len(keys)


def selftest():
    idx = {t: i for i, t in enumerate("ABCDE")}
    # identical topologies -> ASDSF 0
    a = "((A,B),(C,D),E)"
    s1 = newick_splits(a, idx)
    assert frozenset({idx["A"], idx["B"]}) in s1 or frozenset(
        set(range(5)) - {idx["A"], idx["B"]}) in s1, s1
    f = [{sp: 1.0 for sp in s1}] * 3
    m, mx, k = asdsf(f, 0.0)
    assert abs(m) < 1e-12, m
    # one run disagreeing -> ASDSF > 0
    b = "((A,C),(B,D),E)"
    s2 = newick_splits(b, idx)
    m2, _, _ = asdsf([{sp: 1.0 for sp in s1}, {sp: 1.0 for sp in s1},
                      {sp: 1.0 for sp in s2}], 0.0)
    assert m2 > 0.1, m2
    # sampled ancestors (degree-2 node) must not invent a split
    assert newick_splits("((A,B),(C,D),E)", idx) == newick_splits("(((A,B)),(C,D),E)", idx)
    print("selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--burnin", type=float, default=0.10)
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--minfreq", type=float, default=0.10,
                    help="MrBayes-style: splits below this in every run are ignored")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        selftest(); return
    if len(a.files) < 2:
        raise SystemExit("need >= 2 run files (or --selftest)")

    taxa = read_taxa(a.files[0])
    print(f"{len(taxa)} taxa; burnin {a.burnin:.0%}, thin {a.thin}\n")
    freqs = []
    for p in a.files:
        f, n = run_freqs(p, taxa, a.burnin, a.thin)
        freqs.append(f)
        print(f"  {p.split('/')[-1]}: {n} trees, {len(f)} distinct splits")

    m1, mx1, k1 = asdsf(freqs, a.minfreq)
    m0, mx0, k0 = asdsf(freqs, 0.0)
    print(f"\nASDSF (splits >= {a.minfreq:g} in >=1 run, n={k1}) = {m1:.4f}   max SDSF {mx1:.4f}")
    print(f"ASDSF (all observed splits,        n={k0}) = {m0:.4f}   max SDSF {mx0:.4f}")
    v = m1
    verdict = ("EXCELLENT - runs agree on topology" if v < 0.01 else
               "ACCEPTABLE" if v < 0.05 else
               "MARGINAL - treat topology with caution" if v < 0.10 else
               "NOT CONVERGED - the runs are exploring different tree space")
    print(f"\n-> {verdict}")


if __name__ == "__main__":
    main()
