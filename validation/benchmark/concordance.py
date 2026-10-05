#!/usr/bin/env python3
"""
concordance.py — cross-tool CONCORDANCE analysis (no ground truth) for the
GTDB-scale sweep. Measures where ncycle / KofamScan / METABOLIC / NCycDB AGREE
and where they systematically DIVERGE, with a spotlight on the homology-trap
loci (nxrA<->narG, amoA clades, nosZ). DRAM is excluded at scale (see prereg).

We deliberately do NOT score accuracy here: arbitrary GTDB genomes have no
curated truth. Instead we quantify (a) baseline pairwise agreement, (b) the
trap loci where tools split, and (c) ncycle's distinctive calls.

Reuses the validated loaders in adapters.py by pointing them at the pilot
results dir (monkeypatch of module paths) — identical call semantics to the
curated-panel benchmark.

Usage:
  python validation/benchmark/concordance.py \
      --results   comparators/gtdb_pilot/results \
      --ncycdb    comparators/gtdb_pilot/ncycdb.tsv \
      --metabolic comparators/gtdb_pilot/metabolic.tsv \
      --panel-map comparators/gtdb_pilot/panel_map.tsv \
      --out       comparators/gtdb_pilot/CONCORDANCE.md
"""
from __future__ import annotations
import argparse, csv, itertools, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import adapters as A

TOOLS = ["ncycle", "kofam", "metabolic", "ncycdb"]
TRAP_SPOTLIGHT = ["nxrA", "nxrB", "narG", "narH", "amoA", "amoA_archaeal",
                  "amoA_gamma", "nosZ", "nirK", "nrfH"]


def representable() -> dict[str, set[str]]:
    """Targets each tool can in principle call (else a 'disagreement' is just a
    coverage gap, which we report separately)."""
    tkos = A.target_kos()
    has_ko = {t for t, kos in tkos.items() if kos}
    allt = set(A.all_target_ids())
    return {
        "ncycle": allt,
        "kofam": has_ko,
        "metabolic": has_ko,
        "ncycdb": set(A.NCYC_MAP.values()),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--ncycdb", type=Path, help="NCycDB normalized TSV (optional)")
    ap.add_argument("--metabolic", type=Path, help="METABOLIC normalized TSV (optional)")
    ap.add_argument("--panel-map", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    # point the validated loaders at the pilot results dir
    A.RESULTS = args.results
    A.MATRIX = args.results / "ncycle_matrix.tsv"

    genomes = set()
    with open(A.MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genomes.add(r["sample"])

    # ncycle + kofam always available (from the ncycle run); ncycdb/metabolic
    # join only when their normalized TSV exists (partial runs are fine).
    preds = {
        "ncycle": A.load_ncycle(genomes),
        "kofam": A.load_kofam(genomes),
    }
    if args.metabolic and args.metabolic.exists():
        preds["metabolic"] = A.load_metabolic(args.metabolic, genomes)
    if args.ncycdb and args.ncycdb.exists():
        preds["ncycdb"] = A.load_ncycdb(args.ncycdb, genomes)
    global TOOLS
    TOOLS = [t for t in ["ncycle", "kofam", "metabolic", "ncycdb"] if t in preds]
    print(f"[concordance] tools available: {TOOLS}", file=sys.stderr)
    targets = A.all_target_ids()
    repr_ = representable()
    n_g = len(genomes)

    # ── 1. pairwise agreement ────────────────────────────────────────────
    def agreement(a, b, target_set):
        """fraction of (genome,target) cells where tool a and b make the same
        present/absent call, over targets both can represent."""
        tset = [t for t in target_set if t in repr_[a] and t in repr_[b]]
        same = tot = 0
        for g in genomes:
            for t in tset:
                tot += 1
                if preds[a].get((g, t), False) == preds[b].get((g, t), False):
                    same += 1
        return (same / tot) if tot else float("nan"), tot

    def jaccard(a, b, target_set):
        """Positive-call agreement: |present in both| / |present in either|, over
        targets both can represent. Ignores the many shared-absent cells that
        inflate raw agreement — the honest concordance metric."""
        tset = [t for t in target_set if t in repr_[a] and t in repr_[b]]
        both = either = 0
        for g in genomes:
            for t in tset:
                pa = preds[a].get((g, t), False)
                pb = preds[b].get((g, t), False)
                if pa or pb:
                    either += 1
                    if pa and pb:
                        both += 1
        return (both / either) if either else float("nan")

    trap = sorted(A.TRAP)
    pairs = list(itertools.combinations(TOOLS, 2))
    agree_all = {p: agreement(*p, targets) for p in pairs}
    agree_trap = {p: agreement(*p, trap) for p in pairs}
    jac_all = {p: jaccard(*p, targets) for p in pairs}
    jac_trap = {p: jaccard(*p, trap) for p in pairs}

    # ── 2. per-target present counts per tool + divergence ───────────────
    per_target = []
    for t in targets:
        row = {"target": t}
        counts = {}
        for tool in TOOLS:
            if t in repr_[tool]:
                counts[tool] = sum(1 for g in genomes if preds[tool].get((g, t), False))
                row[tool] = counts[tool]
            else:
                row[tool] = "-"   # tool cannot represent this target
        # divergence = spread between max and min present-count among representing tools
        row["divergence"] = (max(counts.values()) - min(counts.values())) if counts else 0
        per_target.append(row)

    # ── 3. ncycle distinctiveness: cells where ncycle differs from ALL others
    distinct = defaultdict(int)   # target -> n genomes where ncycle != every comparator
    for g in genomes:
        for t in targets:
            nc = preds["ncycle"].get((g, t), False)
            others = [preds[o].get((g, t), False) for o in TOOLS[1:] if t in repr_[o]]
            if others and all(o != nc for o in others):
                distinct[t] += 1

    # ── 4. nxrA<->narG trap pattern (the headline locus) ─────────────────
    nxrA_narG = {"ncycle_nxrA_only": 0, "comparator_narG_for_nxr": 0, "rows": []}
    cmp_tools = [t for t in ("kofam", "metabolic", "ncycdb") if t in preds]
    for g in genomes:
        nc_nxrA = preds["ncycle"].get((g, "nxrA"), False)
        nc_narG = preds["ncycle"].get((g, "narG"), False)
        if nc_nxrA:
            nxrA_narG["ncycle_nxrA_only"] += 1
            # which comparators call narG (but not nxrA) on this NOB-like genome?
            for o in cmp_tools:
                if "nxrA" in repr_[o]:
                    if preds[o].get((g, "narG"), False) and not preds[o].get((g, "nxrA"), False):
                        nxrA_narG["comparator_narG_for_nxr"] += 1
                        break

    # ── write report ─────────────────────────────────────────────────────
    labels = {}
    if args.panel_map and args.panel_map.exists():
        for r in csv.DictReader(open(args.panel_map), delimiter="\t"):
            labels[r["sample"]] = r
    n_enriched = sum(1 for r in labels.values() if r.get("stratum") == "enriched")

    L = []
    L.append("# ncycle-pipeline — GTDB-scale Cross-Tool CONCORDANCE (pilot)\n")
    L.append(f"**Genomes:** {n_g} GTDB-representative ({n_enriched} N-cycle-clade-enriched + "
             f"{n_g - n_enriched} cross-phylum backbone) · **Tools:** ncycle, raw KofamScan, "
             "METABOLIC v4.0, NCycDB · **No ground truth** — agreement/divergence only.\n")
    L.append("> Concordance, not accuracy: arbitrary GTDB genomes have no curated truth. "
             "High agreement on specific markers + systematic divergence at the homology traps "
             "is the expected signature; the trap loci are where ncycle's clade resolution acts.\n")

    L.append("\n## 1. Pairwise agreement\n")
    L.append("Over targets both tools can represent. **Raw** = identical present/absent call "
             "(inflated by shared-absent cells). **Jaccard** = positive-call agreement "
             "|both present| / |either present| (the honest concordance metric). Reported for "
             "ALL targets and the homology-trap subset.\n")
    L.append("| tool pair | raw (ALL) | raw (trap) | Jaccard (ALL) | Jaccard (trap) |")
    L.append("|---|---|---|---|---|")
    for p in pairs:
        a_all, _ = agree_all[p]
        a_tr, _ = agree_trap[p]
        L.append(f"| {p[0]} vs {p[1]} | {a_all:.3f} | {a_tr:.3f} | "
                 f"{jac_all[p]:.3f} | {jac_trap[p]:.3f} |")

    L.append("\n## 2. nxrA↔narG trap locus (headline)\n")
    L.append(f"- ncycle calls **nxrA** in **{nxrA_narG['ncycle_nxrA_only']}** genomes.")
    L.append(f"- Of those, **{nxrA_narG['comparator_narG_for_nxr']}** have ≥1 comparator calling "
             "**narG (and NOT nxrA)** on the same genome — the mis-routing of nitrite-oxidizer "
             "NxrA into the nitrate-reductase family that the curated-panel benchmark documented, "
             "now observed at scale.\n")

    L.append("\n## 3. Targets where ncycle is DISTINCTIVE (differs from every comparator)\n")
    L.append("| target | n genomes ncycle-distinct | trap? |")
    L.append("|---|---|---|")
    for t, n in sorted(distinct.items(), key=lambda x: -x[1])[:15]:
        if n:
            L.append(f"| {t} | {n} | {'trap' if t in A.TRAP else ''} |")

    L.append("\n## 4. Per-target present-call counts (by tool; '-' = tool cannot represent)\n")
    L.append("Sorted by divergence (max−min present-count among representing tools).\n")
    L.append("| target | " + " | ".join(TOOLS) + " | divergence |")
    L.append("|---|" + "|".join("---" for _ in TOOLS) + "|---|")
    for row in sorted(per_target, key=lambda r: -r["divergence"]):
        if row["divergence"] == 0 and all(row[t] in (0, "-") for t in TOOLS):
            continue
        L.append("| " + row["target"] + " | "
                 + " | ".join(str(row[t]) for t in TOOLS)
                 + f" | {row['divergence']} |")

    args.out.write_text("\n".join(L) + "\n")
    # also a machine-readable per-target TSV
    tsv = args.out.with_suffix(".tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "ncycle", "kofam", "metabolic",
                                           "ncycdb", "divergence"], delimiter="\t")
        w.writeheader()
        for row in per_target:
            w.writerow(row)
    print(f"[concordance] {n_g} genomes, {len(targets)} targets")
    print(f"[concordance] wrote {args.out} + {tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
