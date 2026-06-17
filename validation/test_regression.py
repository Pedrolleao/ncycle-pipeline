#!/usr/bin/env python3
"""
test_regression.py — accuracy regression gate for ncycle-pipeline.

Fails (exit 1) if panel accuracy drops below locked-in floors. Run AFTER the
pipeline has scored the 33-genome reference panel into
`results/multisample_matrix.tsv` (`make regression` runs pipeline → score → gate).

Floors sit just below the validated curated_v4 baseline (ALL micro-F1 0.97,
precision 0.99, FP 2; homology-trap precision 1.00; hold-out 1.00; min pathway
0.75 = the inherent proteome-mode nitrite-ox limit). They are tight enough to
catch a real regression (a broken HMM, a bad seed edit, a mangled threshold) and
loose enough to tolerate minor scoring noise. Tighten as the tool matures.

Usable two ways:
  python validation/test_regression.py     # prints PASS/FAIL table, exits 0/1
  pytest validation/test_regression.py      # one test per check
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_ncycle import compute_metrics  # noqa: E402

# ── CI lower-bound floors (ROADMAP P5.4.4, 2026-05-31) ───────────────────────
#
# Switched from point-estimate floors to genome-cluster bootstrap CI lower
# bounds (P5.1.1 wired the CIs into score_ncycle.py; the audit's acceptance
# criterion #2). The principle: a real regression is when the CI lower
# bound drops below the established v12 baseline — not when the point
# estimate fluctuates within its already-quantified uncertainty.
#
# Baseline pinned at v12 (post-P5.3.5), with explicit slack below each
# observed CI lower bound to admit normal panel evolution while catching
# a real regression:
#
#   metric              v12 observed CI    floor   slack
#   ALL F1              [1.00, 1.00]       1.00     0     (training is perfect)
#   ALL precision       [1.00, 1.00]       1.00     0     (no training FPs)
#   trap precision      [1.00, 1.00]       0.95     0.05  (traps must stay clean)
#   hold-out F1         [0.695, 0.948]     0.65     0.045 (honest generalization)
#   per-pathway F1      [1.00, 1.00] all   0.70     0.30  (pathway-level safety)
#
# When the CI lower bound is NaN (e.g., empty scope), the check fails.
# MIN_SCORED_CELLS + MAX_ALL_FP stay as point checks (orthogonal to CI).
MIN_ALL_F1_CI_LO          = 1.00
MIN_ALL_PRECISION_CI_LO   = 1.00
MIN_TRAP_PRECISION_CI_LO  = 0.95
MIN_INDEP_TRAP_PRECISION_CI_LO = 0.95   # independent-only (non-circular); observed [1.00,1.00]
MIN_HOLDOUT_F1_CI_LO      = 0.65
MIN_PATHWAY_F1_CI_LO      = 0.70

MAX_ALL_FP                = 4     # orthogonal to CI: an FP regression that
                                  # leaves F1 CI intact is still a regression.
MIN_SCORED_CELLS          = 600   # current 645 training + ~194 hold-out;
                                  # 600 floor guards against an empty/partial run.


def _ci_lo(metric: dict, key: str) -> float:
    """Return CI lower bound for a metric, NaN if missing."""
    import math
    ci = metric.get(f"{key}_ci")
    return ci[0] if ci is not None else float("nan")


def _ci_str(metric: dict, key: str) -> str:
    """Return human-readable point [lo, hi] string."""
    import math
    p = metric[key]
    ci = metric.get(f"{key}_ci")
    if ci is None or math.isnan(ci[0]) or math.isnan(ci[1]):
        return f"{p:.3f}"
    return f"{p:.3f} [{ci[0]:.2f}, {ci[1]:.2f}]"


def _checks(M: dict) -> list[tuple[str, str, bool]]:
    """Return [(check_name, observed_value_str, passed)]."""
    a = M["aggregates"]["ALL"]
    trap = M["aggregates"]["trap"]
    itrap = M["indep_trap"]   # independent-only trap, seed-sourced removed (Frame B)
    # Guard on the DE-LEAKED hold-out (stats audit 2026-06-10): seed-contaminated
    # cells are excluded so the floor tracks the honest generalization number.
    hold = M.get("holdout_deleaked", M["holdout"])
    n_cells = a["TP"] + a["FP"] + a["FN"] + a["TN"]
    checks = [
        (f"scored cells >= {MIN_SCORED_CELLS}",
         str(n_cells), n_cells >= MIN_SCORED_CELLS),
        (f"ALL micro-F1 CI lo >= {MIN_ALL_F1_CI_LO}",
         _ci_str(a, "f1"), _ci_lo(a, "f1") >= MIN_ALL_F1_CI_LO),
        (f"ALL precision CI lo >= {MIN_ALL_PRECISION_CI_LO}",
         _ci_str(a, "precision"), _ci_lo(a, "precision") >= MIN_ALL_PRECISION_CI_LO),
        (f"ALL false-positives <= {MAX_ALL_FP}",
         str(a["FP"]), a["FP"] <= MAX_ALL_FP),
        (f"homology-trap precision CI lo >= {MIN_TRAP_PRECISION_CI_LO}",
         _ci_str(trap, "precision"), _ci_lo(trap, "precision") >= MIN_TRAP_PRECISION_CI_LO),
        (f"trap-independence precision CI lo >= {MIN_INDEP_TRAP_PRECISION_CI_LO}",
         _ci_str(itrap, "precision"), _ci_lo(itrap, "precision") >= MIN_INDEP_TRAP_PRECISION_CI_LO),
        (f"hold-out F1 CI lo >= {MIN_HOLDOUT_F1_CI_LO}",
         _ci_str(hold, "f1"), _ci_lo(hold, "f1") >= MIN_HOLDOUT_F1_CI_LO),
    ]
    for cat, pp in sorted(M["per_pathway"].items()):
        checks.append((f"pathway {cat} F1 CI lo >= {MIN_PATHWAY_F1_CI_LO}",
                       _ci_str(pp, "f1"), _ci_lo(pp, "f1") >= MIN_PATHWAY_F1_CI_LO))
    return checks


# ── pytest entry points (one assertion per check) ────────────────────────────

def test_regression():
    failures = [name for name, _, ok in _checks(compute_metrics()) if not ok]
    assert not failures, "accuracy regression: " + "; ".join(failures)


# ── CLI entry point (human-readable table + exit code) ───────────────────────

def main() -> int:
    results = _checks(compute_metrics())
    width = max(len(name) for name, _, _ in results)
    print("ncycle-pipeline accuracy regression gate (panel = results/multisample_matrix.tsv)\n")
    n_fail = 0
    for name, value, ok in results:
        flag = "PASS" if ok else "FAIL"
        if not ok:
            n_fail += 1
        print(f"  [{flag}] {name:<{width}}  observed={value}")
    print()
    if n_fail:
        print(f"REGRESSION: {n_fail} check(s) failed — accuracy dropped below floor.")
        return 1
    print(f"OK: all {len(results)} checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
