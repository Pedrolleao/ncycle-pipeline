#!/usr/bin/env python3
"""
decompose_gt.py — GT vs detection decomposition (ROADMAP P5.1.3).

Audit ask (REPORT § Audit 2026-05-30, computational-biology review):
separate the GT-revision contribution from the detection-improvement
contribution to the headline F1 trajectory. Full 2×2 factorial would need
both legs:

  (a) hold detection constant, vary GT version:   COMPUTABLE
  (b) hold GT constant, vary detection version:   NOT computable retroactively
      (baseline pipeline state — code + KOfam DB + custom HMMs + BLAST seeds —
      isn't preserved in the repo for each prior point in time). The
      chronological F1 trajectory recorded in REPORT.md historical sections
      is the honest record for this leg.

This script computes leg (a): reconstruct each historical GT (v3 through
v11) by applying the documented cell-flips in reverse, then re-score the
*current* call matrix (results/multisample_matrix.tsv) against each
reconstructed GT. The F1 trajectory across GT versions isolates the
GT-revision contribution.

It also tags every cell-flip with one of three attribution categories so
reviewers can see which moves were biology curation vs detection-enabled:

  - `biology_correction`  — original GT was a curation error; the flip
                            reflects literature / annotation / operon-synteny
                            evidence independent of detection improvements.
  - `detection_enabled`   — original GT was an undetectable-but-real call;
                            the flip became confirmable only after a new HMM,
                            BLAST seed set, or KO fallback added detection
                            machinery. Without that detection improvement
                            the GT would have stayed at the wrong call.
  - `audit_retraction`    — v11 retraction of a prior flip whose original
                            evidence was contaminated (Burkholderia panel
                            file). Cell returned to its v3/v8 state.

Run standalone: prints a markdown-flavored decomposition table that gets
copied into REPORT.md § Audit 2026-05-30 → "Post-audit response (P5.1.3)".
"""
from __future__ import annotations

import math
import sys
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation"))
from build_ground_truth import CURATED, HOLDOUTS, UREASE  # noqa: E402
from score_ncycle import load_calls, prf  # noqa: E402

# ── per-version reverse-transform recipes ────────────────────────────────────
#
# Each FLIPS list captures EVERY documented cell-level change between two
# consecutive GT versions (from build_ground_truth.py docstring + per-genome
# inline comments). To reconstruct an older GT from v11, walk backward
# through the FLIPS lists and apply each in reverse.
#
# Each entry: (genome, target, v_prev, v_next, attribution, note)
# where v_prev / v_next ∈ {"P", "A", None} (None = cell does not exist yet).

FLIPS_V3_TO_V4 = [
    ("Nwinogradskyi_Nb255", "narG", "A", "P", "biology_correction",
     "Starkenburg 2006/2008 + NCBI GFF + operon synteny — v3 was a curation error"),
]
FLIPS_V4_TO_V5 = [
    ("Ngracilis_3211", "nxrB", "P", "A", "biology_correction",
     "GCF_000341545.2 assembly truncation — gene actually absent from this proteome"),
]
FLIPS_V5_TO_V6 = [
    ("Ninopinata_comammox", "nrfA", "A", "P", "biology_correction",
     "CUQ65653.1 annotated by EBI; v5 curation error"),
    ("Nwinogradskyi_Nb255", "nosZ", "A", "P", "audit_retraction",
     "v6 added P based on WP_080511513.1; that accession was Burkholderia "
     "thailandensis NosZ surfacing under the mislabeled panel file (audit "
     "2026-05-30) — RETRACTED in v11 back to A"),
]
FLIPS_V6_TO_V7 = [
    ("Noceani_ATCC19707", "amoA", "P", "A", "detection_enabled",
     "β-AOB amoA HMM genuinely doesn't fire on γ-AOB clade; the new "
     "amoA_gamma HMM made the clade-level call possible — v6 was 'right' "
     "biologically but unscoreable until detection added a γ-clade target"),
    ("Noceani_ATCC19707", "amoA_gamma", None, "P", "detection_enabled",
     "new target added in v7 alongside new amoA_gamma γ-AOB HMM"),
    ("Mfumariolicum_SolV", "amoA_gamma", None, "A", "detection_enabled",
     "primary calibration TN for the new HMM (β-pmoA must not cross-react)"),
    ("Nmultiformis_ATCC25196", "amoA_gamma", None, "A", "detection_enabled",
     "new target TN cell"),
    ("Neuropaea_ATCC19718", "amoA_gamma", None, "A", "detection_enabled",
     "new target TN cell"),
    ("Ninopinata_comammox", "amoA_gamma", None, "A", "detection_enabled",
     "new target TN cell"),
    ("Nnitrosa_comammox", "amoA_gamma", None, "A", "detection_enabled",
     "new target TN cell"),
]
FLIPS_V7_TO_V8 = [
    ("Cnecator_H16", "norB", "P", "A", "biology_correction",
     "Q0JYR9 is qNor (NorB2 locus, 762 aa) not cNorB; original entry comment "
     "already said 'norB qNOR, NO norC' — v3-v7 had a label inconsistency, "
     "v8 closed it"),
]
FLIPS_V8_TO_V9 = [
    ("Cnecator_H16", "norZ", "A", "P", "detection_enabled",
     "Batch 3c qNor BLAST seeds enabled detection of Q7WX97 (Cnecator's "
     "second qNor)"),
    ("Synechocystis_PCC6803", "norZ", "A", "P", "detection_enabled",
     "Batch 3c qNor BLAST seeds enabled detection of P74677 (cyanobacterial "
     "qNor, documented in literature)"),
    ("Nwinogradskyi_Nb255", "norZ", "A", "P", "audit_retraction",
     "v9 added P based on WP_009890104.1 (Burkholderia qNor under mislabeled "
     "panel file); audit 2026-05-30 RETRACTED in v11 — real Nb-255 K04748 "
     "hit BLAST-disqualifies against curated qNor seeds"),
]
FLIPS_V9_TO_V10 = [
    ("Wsuccinogenes_DSM1740", "nosZ", "P", "A", "biology_correction",
     "annotation-gap diagnosis: published clade-II/atypical NosZ carrier "
     "(Simon et al. 2004) but the RefSeq proteome lacks the catalytic NosZ "
     "subunit (zero K00376, zero PF18764, only accessory NosL/NosD present)"),
    # nosZ_clade2 cells added as new target (5 clade-I TN cells); we encode
    # one representative below + a note that the others follow the same pattern.
    ("Pdenitrificans_PD1222", "nosZ_clade2", None, "A", "detection_enabled",
     "new target added in v10; canonical clade-I NosZ TN"),
    ("Paeruginosa_PAO1", "nosZ_clade2", None, "A", "detection_enabled",
     "new target TN cell (γ-proteo clade-I)"),
    ("Cnecator_H16", "nosZ_clade2", None, "A", "detection_enabled",
     "new target TN cell (β-proteo clade-I)"),
    ("Wsuccinogenes_DSM1740", "nosZ_clade2", None, "A", "detection_enabled",
     "new target TN cell; new HMM confirms catalytic NosZ truly absent — "
     "supports the same-version annotation-gap GT flip"),
]
FLIPS_V10_TO_V11 = [
    ("Nwinogradskyi_Nb255", "nosZ", "P", "A", "audit_retraction",
     "P5.0: v6 evidence was Burkholderia; real Nb-255 has only sub-threshold "
     "K00376 (score 36.4) and no PF18764"),
    ("Nwinogradskyi_Nb255", "norZ", "P", "A", "audit_retraction",
     "P5.0: v9 evidence was Burkholderia qNor; real Nb-255 K04748 hit "
     "BLAST-disqualifies against curated qNor seeds"),
    ("Nwinogradskyi_Nb255", "nosZ_clade2", None, "A", "audit_retraction",
     "parallel TN cell added in v11 to keep TN coverage consistent with the "
     "nosZ retraction"),
]

VERSION_FLIPS = [
    ("v3", "v4",  FLIPS_V3_TO_V4),
    ("v4", "v5",  FLIPS_V4_TO_V5),
    ("v5", "v6",  FLIPS_V5_TO_V6),
    ("v6", "v7",  FLIPS_V6_TO_V7),
    ("v7", "v8",  FLIPS_V7_TO_V8),
    ("v8", "v9",  FLIPS_V8_TO_V9),
    ("v9", "v10", FLIPS_V9_TO_V10),
    ("v10","v11", FLIPS_V10_TO_V11),
]
ALL_FLIPS = [(prev, nxt, f) for prev, nxt, fl in VERSION_FLIPS for f in [fl]]


# ── reconstruct a historical GT by applying reverse flips ────────────────────

def build_v11_cells() -> dict[tuple[str, str], str]:
    """v11 ground-truth cells as a flat (genome,target) → 'present'|'absent' dict.
    Replicates the dict that build_ground_truth.py writes to ground_truth.tsv."""
    out = {}
    for genome, sets in CURATED.items():
        for t in dict.fromkeys(sets.get("P", [])):
            out[(genome, t)] = "present"
        for t in dict.fromkeys(sets.get("A", [])):
            out[(genome, t)] = "absent"
    return out


def reverse_flip(cells: dict, flip) -> dict:
    """Apply one cell-flip in reverse (newer-version cell → older-version cell)."""
    out = deepcopy(cells)
    g, t, v_prev, v_next, _attr, _note = flip
    key = (g, t)
    if v_prev is None:
        # Cell didn't exist in the prior version → drop it.
        out.pop(key, None)
    else:
        out[key] = {"P": "present", "A": "absent"}[v_prev]
    return out


def gt_at_version(version: str) -> dict:
    """Build the GT cells dict at any tagged version v3..v11 by walking
    backward from v11 through the reverse flips."""
    if version == "v11":
        return build_v11_cells()
    cells = build_v11_cells()
    for prev, nxt, flips in reversed(VERSION_FLIPS):
        for fl in flips:
            cells = reverse_flip(cells, fl)
        if prev == version:
            return cells
    raise ValueError(f"unknown version: {version}")


# ── score current matrix against any GT ──────────────────────────────────────

def score_against(gt_cells: dict) -> dict:
    """Score the current results/multisample_matrix.tsv predictions against
    a reconstructed GT. Reports training + hold-out aggregates."""
    calls = load_calls()
    train = {"TP":0,"FP":0,"FN":0,"TN":0}
    hold  = {"TP":0,"FP":0,"FN":0,"TN":0}
    for (g, t), exp in gt_cells.items():
        pred_present = calls.get((g, t), 0) in (1, 2)
        exp_pos = (exp == "present")
        cell = ("TP" if pred_present and exp_pos else
                "FP" if pred_present and not exp_pos else
                "FN" if not pred_present and exp_pos else "TN")
        (hold if g in HOLDOUTS else train)[cell] += 1
    tp, rp, tf = prf(train["TP"], train["FP"], train["FN"])
    hp, hr, hf = prf(hold["TP"],  hold["FP"],  hold["FN"])
    return {
        "train": {**train, "precision": tp, "recall": rp, "f1": tf},
        "hold":  {**hold,  "precision": hp, "recall": hr, "f1": hf},
        "n_cells": len(gt_cells),
    }


# ── print: trajectory table + attribution table ──────────────────────────────

def fmt(v): return "NA" if (isinstance(v, float) and math.isnan(v)) else f"{v:.3f}"


def print_trajectory_table():
    print("\n## GT-revision trajectory (current detection held constant)\n")
    print(f"{'GT version':<12}{'n_cells':>10}{'TP':>6}{'FP':>6}{'FN':>6}{'TN':>6}"
          f"{'P':>9}{'R':>9}{'F1':>9}  hold-out F1")
    print("-" * 92)
    for v in ["v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11"]:
        gt = gt_at_version(v)
        s = score_against(gt)
        tr = s["train"]; ho = s["hold"]
        print(f"{v:<12}{s['n_cells']:>10}{tr['TP']:>6}{tr['FP']:>6}{tr['FN']:>6}"
              f"{tr['TN']:>6}{fmt(tr['precision']):>9}{fmt(tr['recall']):>9}"
              f"{fmt(tr['f1']):>9}  {fmt(ho['f1'])}")


def print_attribution_table():
    print("\n## Cell-flip attribution (v3 → v11)\n")
    print(f"{'from':<5}{'to':<5}  {'genome':<28}{'target':<14}{'change':<8}"
          f"{'attribution':<24}note")
    print("-" * 140)
    for prev, nxt, flips in VERSION_FLIPS:
        for (g, t, v_prev, v_next, attr, note) in flips:
            if v_prev is None: change = f"new={v_next}"
            elif v_next is None: change = f"drop"
            else: change = f"{v_prev}→{v_next}"
            print(f"{prev:<5}{nxt:<5}  {g:<28}{t:<14}{change:<8}{attr:<24}{note[:60]}")


def print_summary():
    counts = {"biology_correction": 0, "detection_enabled": 0, "audit_retraction": 0}
    for prev, nxt, flips in VERSION_FLIPS:
        for (_g, _t, _vp, _vn, attr, _note) in flips:
            counts[attr] = counts.get(attr, 0) + 1
    total = sum(counts.values())
    print("\n## Attribution summary (v3 → v11)\n")
    for k, n in counts.items():
        print(f"  {k:<24}{n:>3}/{total} cells")


def main() -> int:
    print("# GT vs detection decomposition  (ROADMAP P5.1.3)")
    print("# Current detection held constant; GT version varies. The leg")
    print("# 'GT held constant, detection varies' is NOT computable retroactively")
    print("# without preserved baseline pipeline state — see REPORT.md")
    print("# historical sections (P1, P2, P3, P3b, Panel expansion, P4, Batches")
    print("# 1/2/3) for the chronological detection-improvement trajectory.")
    print_trajectory_table()
    print_summary()
    print_attribution_table()
    return 0


if __name__ == "__main__":
    sys.exit(main())
