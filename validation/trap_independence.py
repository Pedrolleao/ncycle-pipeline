#!/usr/bin/env python3
"""
trap_independence.py — audit the homology-trap claim for seed/HMM leakage (Frame B).

Ported from scycle-pipeline/validation/trap_independence.py to harmonize the two pipelines'
reporting (both now report: whole-panel F1, an independent-subset metric, and a train/hold-out
split + LOGO-CV). The benchmark's confirmatory endpoint is trap PRECISION. Some training-panel
genomes are also curated-BLAST-seed or custom-HMM sources for the very trap targets they test,
so their positive calls are partly circular (the model was built to detect those proteins).

Independence is a per-(genome, target) property, NOT per-genome: e.g. Escherichia coli K-12 is a
narG/narH/napA/nirB seed source (non-independent there) but an independent test positive for any
trap target it was not used to build. Critically, trap PRECISION is robust to seed attribution —
false positives fall on GT-absent cells, which are never "seeded positives" — so the
INDEPENDENT-only trap precision below is the load-bearing, non-circular number.

Run from the pipeline root:
    python validation/trap_independence.py

SEED_SRC is derived from the actual positive seed/HMM-ref organisms, mapped to TRAINING-panel
genome names:
  - targets/<t>/manifest.yaml kept `expanded_candidates` (custom-HMM training refs), and
  - config/targets.yaml `blast_refs_uniprot`, and resources/blast_db/blast_gated_refs.fasta
    headers (the actual BLAST-gating DB).
Hold-out genomes (source=holdout_v3) are excluded from the training score by load_truth, so only
training-panel attributions matter here.
"""
from __future__ import annotations
import sys
from pathlib import Path

# NB: the benchmark `adapters` import is deferred into main() so that
# `from trap_independence import is_seed, SEED_SRC` is cheap and side-effect-free
# (the regression gate imports those two for its trap-independence row).

# trap target -> TRAINING-panel genomes that are POSITIVE curated-seed / custom-HMM sources for it.
# Derivation (2026-06-07): see header. Only positive seeds create positive-circularity; negatives
# and KO-only targets do not. Targets omitted below have NO training-panel seed source -> every
# panel positive is independent (notably nxrA: refs are N. alkalicus / Nitrospira nitrificans /
# kreftii — none on the panel — so winogradskyi/gracilis/inopinata/japonica/nitrosa are all
# independent nxrA test positives; amoC: Nitrosomonas sp. TK794 + uncultured; amoA_archaeal:
# PF12942 archaeal-specific Pfam, no curated panel seed; nirA: KO-based assimilatory).
SEED_SRC: dict[str, set[str]] = {
    "amoA": {"Neuropaea_ATCC19718", "Nmultiformis_ATCC25196", "Nnitrosa_comammox"},
    "amoB": {"Neuropaea_ATCC19718"},                      # gated ref Q04508 = N. europaea ATCC19718
    "nxrB": {"Ngracilis_3211", "Nnitrosa_comammox"},      # HMM refs M1KVL1 (Ngracilis) + Nitrospira nitrosa
    "narG": {"Ecoli_K12_MG1655", "Bsubtilis_168"},        # gated refs P09152 (Eco) + P42175 (Bsu)
    "narH": {"Ecoli_K12_MG1655", "Bsubtilis_168"},        # gated refs P11349 (Eco) + P42176 (Bsu)
    "napA": {"Ecoli_K12_MG1655", "Cnecator_H16"},         # seeds P33937 (Eco) + P39185 (Cnecator H16)
    "nirB": {"Ecoli_K12_MG1655"},                         # gated ref P08201 (Eco); Bsu entry is NASD (assim)
    "nirD": {"Ecoli_K12_MG1655"},                         # gated ref P0A9I8 (Eco)
    # nirK: blast_seeds (P25006/P38501/Q06006/Q53239/A0A5C7VYN2/A0A1I4L3J7) are NOT resolved to any
    #   panel genome (KO-based target, no curated panel-derived seed confirmed); panel nirK positives
    #   are treated as independent pending a seed-organism audit. Precision is unaffected regardless
    #   (FPs land on absent cells). Flagged here for honesty, mirroring scycle's phsA/ttrA caveat.
}


def is_seed(g: str, t: str) -> bool:
    return g in SEED_SRC.get(t, set())


def _metrics(truth, pred, cells):
    tp = fp = fn = 0
    for (g, t) in cells:
        gt = truth[(g, t)] == "present"
        pr = pred.get((g, t), False)
        tp += pr and gt
        fp += pr and not gt
        fn += (not pr) and gt
    P = tp / (tp + fp) if tp + fp else float("nan")
    R = tp / (tp + fn) if tp + fn else float("nan")
    F = 2 * P * R / (P + R) if (P == P and R == R and P + R) else float("nan")
    return tp, fp, fn, P, R, F


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent / "benchmark"))
    import adapters
    truth, genomes = adapters.load_truth()
    pred = adapters.load_ncycle(genomes)
    TRAP = adapters.TRAP

    present = sorted((g, t) for (g, t), v in truth.items() if t in TRAP and v == "present")
    print(f"=== {len(present)} PRESENT trap cells (GT, training panel), tagged ===")
    for g, t in present:
        tag = "SEED-SOURCED" if is_seed(g, t) else "independent"
        pr = "present" if pred.get((g, t), False) else "ABSENT"
        print(f"  {t:14} {g:24} pred={pr:8} [{tag}]")

    trap_cells = [(g, t) for (g, t) in truth if t in TRAP]
    indep_cells = [c for c in trap_cells if not is_seed(*c)]
    print("\n=== ncycle trap metrics ===")
    for label, cells in [("ALL trap", trap_cells), ("INDEPENDENT-only trap", indep_cells)]:
        tp, fp, fn, P, R, F = _metrics(truth, pred, cells)
        npos = sum(truth[c] == "present" for c in cells)
        print(f"  {label:24} cells={len(cells)} present={npos}  "
              f"TP={tp} FP={fp} FN={fn}  P={P:.3f} R={R:.3f} F1={F:.3f}")


if __name__ == "__main__":
    main()
