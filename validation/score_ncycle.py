#!/usr/bin/env python3
"""
score_ncycle.py — score pipeline calls against the curated N-cycle ground truth.

Adapted from ewaste-pipeline/validation/scripts/score_phase1.py. Differences:
  - scores ONLY the (genome,target) cells present in ground_truth.tsv (curated
    subset), instead of defaulting every unlisted cell to absent;
  - reports per-pathway breakdown + a homology-trap subset instead of Tier-1/2;
  - reports 95% genome-cluster bootstrap CIs alongside each aggregate / per-
    pathway / hold-out point estimate (ROADMAP P5.1.1, added 2026-05-30 to
    meet acceptance criterion #2). Method mirrors `validation/benchmark/
    benchmark_stats.py:130-144` `boot_ci`: resample genomes with replacement
    B=10,000 times, recompute the metric on each resample, take 2.5/97.5
    percentiles. The unit of independence is the genome (cells within a
    genome are correlated).

Matrix status codes (ncycle_results/ncycle_matrix.tsv): 2 confirmed, 1 domain-only
or narrow-no-IPR (→ predicted present), 0 absent, -1 disqualified (→ absent).
"""
from __future__ import annotations
import csv, math, os, random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Default GT = manual/literature-curated ground truth (ground_truth.tsv, the
# authoritative reference). Set env GT_FILE=<file> to score against an alternate
# GT in validation/ (uniform interface with scycle-pipeline's GT_FILE switch — used
# for the automated-KEGG-GT contrast once that second GT is built; see WORKPLAN B2).
GT = ROOT / "validation" / os.environ.get("GT_FILE", "ground_truth.tsv")
MATRIX = ROOT / "ncycle_results" / "ncycle_matrix.tsv"
TARGETS = ROOT / "config" / "targets.yaml"
# one pair of output tables per ground truth, so that scoring the KEGG contrast does
# not overwrite the tables of the default (manual) ground truth
_TAG = "" if GT.name == "ground_truth.tsv" else "." + GT.stem
OUT_M = ROOT / "validation" / f"ncycle_metrics{_TAG}.tsv"
OUT_C = ROOT / "validation" / f"ncycle_confusion{_TAG}.tsv"
LEAKAGE = ROOT / "validation" / "holdout_seed_leakage.tsv"

TRAP = {"nxrA","nxrB","narG","narH","amoA","amoB","amoC","napA",
        "nirK","nirB","nirD","nirA","amoA_archaeal"}

# Bootstrap CI parameters (genome-cluster, percentile method).
BOOTSTRAP_B    = 10_000   # B=10k matches benchmark_stats.py
BOOTSTRAP_SEED = 1234     # deterministic across runs

# Per-target F1 denominator floor (ROADMAP P5.1.4). A target whose panel has
# fewer than this many positive-class cells (TP+FN) produces an F1 estimate
# dominated by 1-cell granularity; we flag it as low-reliability in both
# stdout (asterisk + legend) and the TSV (boolean f1_reliable column). The
# point estimate is still computed and reported — just annotated.
MIN_RELIABLE_POSITIVES = 3


def load_truth():
    out = {}
    holdouts = set()
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["genome"], r["target"])] = r["expected"]
            if r.get("source") == "holdout_v3":
                holdouts.add(r["genome"])
    return out, holdouts


def load_leakage():
    """Set of (genome, target) hold-out cells excluded as seed-contaminated — the
    target's detector (BLAST gate or custom HMM) was seeded from a conspecific/
    congeneric organism, so the cell measures memorization not generalization.
    Produced by detect_seed_leakage.py; empty set if the file is absent."""
    out = set()
    if LEAKAGE.exists():
        with open(LEAKAGE) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.add((r["genome"], r["target"]))
    return out


def load_calls():
    out = {}
    with open(MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["sample"]
            for col, val in r.items():
                if col.startswith("target__"):
                    try: out[(g, col[len("target__"):])] = int(val)
                    except (TypeError, ValueError): out[(g, col[len("target__"):])] = 0
    return out


def cat_map():
    import yaml
    cfg = yaml.safe_load(open(TARGETS))
    return {t["id"]: t["category"] for t in cfg["targets"]}


def prf(tp, fp, fn):
    p = tp/(tp+fp) if (tp+fp) else float("nan")
    r = tp/(tp+fn) if (tp+fn) else float("nan")
    f = (2*p*r/(p+r)) if (not math.isnan(p) and not math.isnan(r) and (p+r)) else float("nan")
    return p, r, f


def rule_of_three_upper(n_cells: int) -> float:
    """One-sided ~95% upper bound on the per-cell error rate given 0 observed
    errors in n_cells (rule of three: 3/n). Replaces the degenerate [point,point]
    bootstrap CI for a zero-error scope: '0 errors is consistent with a true error
    rate up to X%', NOT 'the error rate is exactly 0'."""
    return 3.0 / n_cells if n_cells else float("nan")


def _sum_tallies(tallies):
    """Sum a list of (TP,FP,FN,TN) tuples element-wise."""
    if not tallies:
        return (0, 0, 0, 0)
    return tuple(sum(x) for x in zip(*tallies))


def _metric_of(tally, name):
    """Compute precision|recall|f1 from a single (TP,FP,FN,TN) tuple."""
    p, r, f = prf(tally[0], tally[1], tally[2])
    return {"precision": p, "recall": r, "f1": f}[name]


def boot_ci(gtallies: dict, name: str, *,
            B: int = BOOTSTRAP_B,
            rng: random.Random | None = None) -> tuple[float, float]:
    """Genome-cluster 95% percentile bootstrap CI for precision|recall|f1.

    `gtallies`: dict[genome -> (TP, FP, FN, TN)] over the scope of interest.
    Resamples genomes with replacement B times, recomputes the metric on each
    resampled set, returns (2.5th, 97.5th) percentiles. NaN-only resamples
    (e.g. all-TN with no positives so recall is undefined) are dropped before
    percentile estimation. Mirrors `validation/benchmark/benchmark_stats.py:
    130-144` so headline-path CIs use the same method as the benchmark suite.

    NOTE: when every panel genome scores cleanly (FP=0 AND FN=0 across the
    scope), every resample also scores cleanly → the CI collapses to
    [point, point]. This reflects that the cluster bootstrap cannot
    manufacture variability from a panel that has none; it does NOT mean the
    estimate is precise in an absolute sense (for a zero-error scope, see
    `rule_of_three_upper` for the honest upper bound on the error rate). The
    hold-out (8 genomes) gives a non-degenerate interval; the training panel,
    being calibration data, typically scores cleanly and collapses.
    """
    if rng is None:
        rng = random.Random(BOOTSTRAP_SEED)
    glist = list(gtallies)
    if not glist:
        return (float("nan"), float("nan"))
    vals = []
    for _ in range(B):
        samp = [gtallies[rng.choice(glist)] for _ in glist]
        v = _metric_of(_sum_tallies(samp), name)
        if not math.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[min(len(vals) - 1, int(0.975 * len(vals)))]
    return (lo, hi)


def compute_metrics() -> dict:
    """Score the current matrix against ground truth and return a structured
    metrics dict (single source of truth for both the stdout report in main()
    and the regression gate in test_regression.py). No file I/O / printing.

    As of P5.1.1 (2026-05-30) each aggregate / per-pathway / hold-out entry
    also carries genome-cluster bootstrap 95% CIs for precision, recall, F1
    (where applicable). See `boot_ci()` for method + the all-clean-collapse
    caveat.
    """
    truth, holdouts = load_truth(); calls = load_calls(); cats = cat_map()
    rng = random.Random(BOOTSTRAP_SEED)

    # First pass: per-cell records + confusion table (training vs hold-out split).
    train_cells: list[tuple] = []   # (genome, target, category, cell_label)
    hold_cells: list[tuple] = []
    confusion = []
    for (g, t), exp in truth.items():
        pred_code = calls.get((g, t), 0)
        pred_present = pred_code in (1, 2)
        exp_pos = (exp == "present")
        cell = ("TP" if pred_present and exp_pos else
                "FP" if pred_present and not exp_pos else
                "FN" if not pred_present and exp_pos else "TN")
        confusion.append((g, t, exp, pred_code,
                          "present" if pred_present else "absent", cell))
        rec = (g, t, cats.get(t, "?"), cell)
        (hold_cells if g in holdouts else train_cells).append(rec)

    # Per-target tallies (training only) — drives per-target output rows.
    per: dict[str, dict[str, int]] = {}
    for (g, t, cat, cell) in train_cells:
        d = per.setdefault(t, {"TP":0,"FP":0,"FN":0,"TN":0})
        d[cell] += 1

    rows = []
    for t, d in per.items():
        p, r, f = prf(d["TP"], d["FP"], d["FN"])
        n_pos = d["TP"] + d["FN"]
        rows.append({"target": t, "category": cats.get(t,"?"),
                     "n": sum(d.values()), "n_pos": n_pos,
                     "f1_reliable": n_pos >= MIN_RELIABLE_POSITIVES,
                     **d,
                     "precision": p, "recall": r, "f1": f})

    # ── per-genome tallies under each scope (cluster-bootstrap unit) ────────
    _IDX = {"TP":0, "FP":1, "FN":2, "TN":3}

    def _gtallies(scope_cells):
        out: dict[str, list[int]] = {}
        for (g, _t, _cat, cell) in scope_cells:
            d = out.setdefault(g, [0, 0, 0, 0])
            d[_IDX[cell]] += 1
        return {g: tuple(v) for g, v in out.items()}

    train_trap_cells    = [c for c in train_cells if c[1] in TRAP]
    train_nontrap_cells = [c for c in train_cells if c[1] not in TRAP]
    pathway_cells: dict[str, list[tuple]] = {}
    for cat in sorted({c[2] for c in train_cells}):
        pathway_cells[cat] = [c for c in train_cells if c[2] == cat]

    # ── aggregate builder (point + CIs) ─────────────────────────────────────
    def aggregate(scope_cells, *, tids_for_median: list[str] | None = None,
                  with_pr_ci: bool = True) -> dict:
        gt = _gtallies(scope_cells)
        TP, FP, FN, TN = _sum_tallies(list(gt.values()))
        p, r, f = prf(TP, FP, FN)
        out = {"TP": TP, "FP": FP, "FN": FN, "TN": TN,
               "precision": p, "recall": r, "f1": f,
               "f1_ci": boot_ci(gt, "f1", rng=rng)}
        if with_pr_ci:
            out["precision_ci"] = boot_ci(gt, "precision", rng=rng)
            out["recall_ci"]    = boot_ci(gt, "recall",    rng=rng)
        if tids_for_median is not None:
            f1s = sorted(m["f1"] for m in rows
                         if m["target"] in tids_for_median and not math.isnan(m["f1"]))
            out["median_f1"] = f1s[len(f1s)//2] if f1s else float("nan")
        return out

    aggregates = {
        "ALL":     aggregate(train_cells,         tids_for_median=list(per)),
        "trap":    aggregate(train_trap_cells,    tids_for_median=[t for t in per if t in TRAP]),
        "nontrap": aggregate(train_nontrap_cells, tids_for_median=[t for t in per if t not in TRAP]),
    }

    # ── Frame B: independent-only trap precision (harmonized gate row, 2026-06-12) ─
    # Trap precision over training-panel trap cells EXCLUDING seed/HMM-sourced
    # (g,t) pairs — the non-circular number. is_seed/SEED_SRC from trap_independence.py.
    from trap_independence import is_seed
    indep_trap = aggregate([c for c in train_trap_cells if not is_seed(c[0], c[1])])

    per_pathway = {}
    for cat in sorted(pathway_cells):
        agg = aggregate(pathway_cells[cat], with_pr_ci=False)
        per_pathway[cat] = {k: agg[k] for k in ("TP", "FP", "FN", "f1", "f1_ci")}

    holdout = aggregate(hold_cells)
    holdout["genomes"] = sorted(holdouts)

    # ── de-leaked hold-out (drop seed-contaminated cells) ───────────────────
    # The headline generalization number: excludes cells whose detector was
    # seeded from a conspecific/congeneric organism (detect_seed_leakage.py),
    # which trivially pass and measure memorization, not generalization.
    leak = load_leakage()
    hold_cells_dl = [c for c in hold_cells if (c[0], c[1]) not in leak]
    holdout_deleaked = aggregate(hold_cells_dl)
    holdout_deleaked["genomes"] = sorted({c[0] for c in hold_cells_dl})
    holdout_deleaked["excluded"] = sorted(leak)

    # ── Frame A: whole-panel (training + hold-out combined) ─────────────────
    # Harmonization (2026-06-07): report whole-panel micro-F1 alongside the
    # train/hold-out split (Frame C) and the independent-subset metric
    # (validation/trap_independence.py, Frame B), mirroring scycle-pipeline.
    # This is the directly-comparable single-number cut; the train/hold-out
    # split above remains the generalization claim, the hold-out CI the guard.
    whole_cells = train_cells + hold_cells
    whole = {
        "ALL":     aggregate(whole_cells),
        "trap":    aggregate([c for c in whole_cells if c[1] in TRAP]),
        "nontrap": aggregate([c for c in whole_cells if c[1] not in TRAP]),
    }
    whole["genomes"] = sorted({c[0] for c in whole_cells})

    return {"per_target": rows, "per": per, "aggregates": aggregates,
            "per_pathway": per_pathway, "indep_trap": indep_trap,
            "holdout": holdout, "holdout_deleaked": holdout_deleaked,
            "whole": whole, "confusion": confusion}


def main():
    M = compute_metrics()
    rows = M["per_target"]; per = M["per"]; confusion = M["confusion"]

    def fmt(v): return "NA" if (isinstance(v,float) and math.isnan(v)) else f"{v:.2f}"

    OUT_M.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_M, "w") as fh:
        fh.write("target\tcategory\tn\tn_pos\tTP\tFP\tFN\tTN\tprecision\trecall\tf1\tf1_reliable\n")
        for m in sorted(rows, key=lambda x: (x["category"], x["target"])):
            fh.write(f"{m['target']}\t{m['category']}\t{m['n']}\t{m['n_pos']}\t"
                     f"{m['TP']}\t{m['FP']}\t{m['FN']}\t{m['TN']}\t"
                     f"{fmt(m['precision'])}\t{fmt(m['recall'])}\t{fmt(m['f1'])}\t"
                     f"{str(m['f1_reliable']).lower()}\n")
    with open(OUT_C, "w") as fh:
        fh.write("genome\ttarget\texpected\tstatus_code\tpredicted\tcell\n")
        for c in confusion:
            fh.write("\t".join(str(x) for x in c) + "\n")

    # ---- stdout report ----
    print(f"{'target':<16}{'pathway':<24}{'TP':>3}{'FP':>3}{'FN':>3}{'TN':>3}  {'P':>4}{'R':>5}{'F1':>6}")
    low_n_targets: list[tuple[str, int]] = []
    for m in sorted(rows, key=lambda x:(x["category"], -(x["TP"]+x["FN"]), x["target"])):
        flag = "*" if not m["f1_reliable"] else " "
        if not m["f1_reliable"]:
            low_n_targets.append((m["target"], m["n_pos"]))
        print(f"  {m['target']:<14}{m['category']:<24}{m['TP']:>3}{m['FP']:>3}{m['FN']:>3}{m['TN']:>3}  "
              f"{fmt(m['precision'])}{fmt(m['recall']):>5}{fmt(m['f1']):>5}{flag}")

    if low_n_targets:
        print(f"\n[*] F1 estimated from <{MIN_RELIABLE_POSITIVES} positive-class cells "
              f"(TP+FN); single FP/FN moves it by 1/n_pos. Interpret with caution. "
              f"Targets affected:")
        for t, n in sorted(low_n_targets, key=lambda x: (x[1], x[0])):
            print(f"      {t} (n_pos={n})")

    def fmt_ci(ci):
        if ci is None: return ""
        lo, hi = ci
        if math.isnan(lo) or math.isnan(hi): return ""
        return f" [{lo:.2f}, {hi:.2f}]"

    def agg(label, a):
        med = (f"  median-target-F1={fmt(a['median_f1'])}" if "median_f1" in a else "")
        n_cells = a['TP'] + a['FP'] + a['FN'] + a['TN']
        # For a zero-error scope the bootstrap CI is degenerate ([1.00,1.00]); report
        # the rule-of-three upper bound on the true error rate instead (stats audit
        # 2026-06-10) so a perfect score is never presented as zero-uncertainty.
        perfect = (a['FP'] == 0 and a['FN'] == 0 and n_cells > 0)
        bound = (f"  [0 err/{n_cells} cells → true per-cell error ≤ "
                 f"{100*rule_of_three_upper(n_cells):.2f}% (rule of 3, 95%); "
                 f"bootstrap CI degenerate]" if perfect else "")
        print(f"\n[{label}] cells={n_cells}  "
              f"TP={a['TP']} FP={a['FP']} FN={a['FN']} TN={a['TN']}  "
              f"micro-P={fmt(a['precision'])}{fmt_ci(a.get('precision_ci'))} "
              f"micro-R={fmt(a['recall'])}{fmt_ci(a.get('recall_ci'))} "
              f"micro-F1={fmt(a['f1'])}{fmt_ci(a.get('f1_ci'))}{med}{bound}")

    print("\n=== Frame C — TRAINING panel (generalization split) ===")
    agg("ALL curated targets", M["aggregates"]["ALL"])
    agg("Homology-trap targets", M["aggregates"]["trap"])
    agg("Non-trap targets", M["aggregates"]["nontrap"])
    print("\nper-pathway micro-F1 (training panel):")
    for cat, pp in sorted(M["per_pathway"].items()):
        print(f"  {cat:<26} F1={fmt(pp['f1'])}{fmt_ci(pp.get('f1_ci'))}  "
              f"(TP={pp['TP']} FP={pp['FP']} FN={pp['FN']})")

    # ---- Frame A: whole-panel (training + hold-out combined) ----
    w = M["whole"]
    print(f"\n=== Frame A — WHOLE PANEL ({len(w['genomes'])} genomes, training + hold-out) ===")
    agg("ALL curated targets", w["ALL"])
    agg("Homology-trap targets", w["trap"])
    agg("Non-trap targets", w["nontrap"])
    print("  (Frame B — independent-subset trap precision: run validation/trap_independence.py)")

    # ---- Frame C: independent HOLD-OUT test set (de-leaked = headline) ----
    hd = M["holdout_deleaked"]; h = M["holdout"]
    excl = hd.get("excluded", [])
    print(f"\n=== Frame C — HOLD-OUT test set ({len(hd['genomes'])} genomes, "
          f"the generalization claim) ===")

    def hline(tag, x):
        print(f"  {tag}: cells={x['TP']+x['FP']+x['FN']+x['TN']}  "
              f"TP={x['TP']} FP={x['FP']} FN={x['FN']} TN={x['TN']}  "
              f"P={fmt(x['precision'])}{fmt_ci(x.get('precision_ci'))} "
              f"R={fmt(x['recall'])}{fmt_ci(x.get('recall_ci'))} "
              f"F1={fmt(x['f1'])}{fmt_ci(x.get('f1_ci'))}")

    hline("DE-LEAKED (headline)     ", hd)
    hline("as-is (seed-contaminated)", h)
    if excl:
        print(f"  excluded {len(excl)} seed-contaminated positive cells "
              f"(detect_seed_leakage.py): " + ", ".join(f"{g}/{t}" for g, t in excl))

    print(f"\n[bootstrap] 95% CIs via genome-cluster percentile bootstrap, "
          f"B={BOOTSTRAP_B}, seed={BOOTSTRAP_SEED}. CIs collapse to [point, point] "
          f"when the scope has zero FP and zero FN (see the rule-of-3 bound printed "
          f"inline for those scopes) — variability can't be manufactured from a panel "
          f"that has none. The training panel is used for threshold calibration, so its "
          f"perfect score is an IN-SAMPLE fit; the de-leaked HOLD-OUT above is the "
          f"credible generalization estimate.")


if __name__ == "__main__":
    main()
