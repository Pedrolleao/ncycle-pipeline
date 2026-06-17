#!/usr/bin/env python3
"""
compute_complex_completeness.py — given calls/ncycle_calls.tsv and the
complexes/synergies blocks of targets.yaml, emit:

  calls/complex_completeness.tsv
      complex_id  completeness  status  members_present  members_missing  application

  calls/synergy_completeness.tsv
      synergy_id  completeness  status  requires_present  requires_missing  benefit
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml

PRESENT_STATUSES = {"confirmed", "domain-only", "narrow-no-IPR"}


def load_calls(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["target_id"]] = row["status"]
    return out


def status_for(
    completeness: float, forbids_violated: int = 0,
    n_slots_filled: int = -1, n_slots_total: int = 0,
    single_organism: bool = False,
) -> str:
    # Hard veto 1: any forbid violation → phenotype categorically does not apply
    # (a complete denitrifier is the OPPOSITE of n2o_emitter, not "partially one").
    if forbids_violated > 0 and completeness < 0.999:
        return "absent"
    # Hard veto 2: if the synergy has required slots and NONE are filled, the
    # phenotype is absent regardless of how many forbids are vacuously satisfied
    # (e.g., a Synechocystis with no N-dissim genes isn't "partial n2o_sink").
    if n_slots_total > 0 and n_slots_filled == 0:
        return "absent"
    if completeness >= 0.999:
        return "complete"
    # Single-organism phenotypes (comammox, anammox) cannot be community-
    # completed — partial is meaningless, so collapse to absent.
    if single_organism:
        return "absent"
    if completeness >= 0.5:
        return "partial"
    return "absent"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--complexes-out", required=True, type=Path)
    ap.add_argument("--synergies-out", required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    complexes = cfg.get("complexes", {})
    synergies = cfg.get("synergies", [])
    calls = load_calls(args.calls)

    # Complexes. Each entry in `members` is a slot: a single target id (string)
    # or a list of equivalents (any-of — paralog substitutes count). Mirrors
    # the synergy evaluator below (ROADMAP P5.3.3, 2026-05-30; ported the
    # any-of slot pattern from the synergy block so complex members can also
    # accept clade paralogs — e.g. ammonia_monooxygenase: [[amoA, amoA_gamma,
    # amoA_archaeal], amoB, amoC] lets γ-AOB Nitrosococcus + AOA complete the
    # complex without per-clade parallel complex entries).
    args.complexes_out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.complexes_out, "w") as fh:
        fh.write("complex_id\tcompleteness\tstatus\tmembers_present\t"
                 "members_missing\tn_present\tn_total\tapplication\n")
        for cid, body in complexes.items():
            slots = [[m] if isinstance(m, str) else list(m)
                     for m in body["members"]]
            present_per_slot = [
                [m for m in slot if calls.get(m) in PRESENT_STATUSES]
                for slot in slots
            ]
            present = [p[0] for p in present_per_slot if p]
            missing = [slot[0] if len(slot) == 1 else "|".join(slot)
                       for slot, p in zip(slots, present_per_slot) if not p]
            n_slots_filled = sum(1 for p in present_per_slot if p)
            comp = n_slots_filled / len(slots) if slots else 0.0
            fh.write("\t".join([
                cid, f"{comp:.3f}", status_for(comp),
                ",".join(present), ",".join(missing),
                str(n_slots_filled), str(len(slots)),
                body.get("application", ""),
            ]) + "\n")

    # Synergies. Each entry in `requires` is a slot: either a single target id
    # (string) or a list of equivalents (any-of — paralog substitutes count).
    # Optional `forbids:` lists targets that must be ABSENT (negative conjunction).
    with open(args.synergies_out, "w") as fh:
        fh.write("synergy_id\tcompleteness\tstatus\trequires_present\t"
                 "requires_missing\tforbids_satisfied\tforbids_violated\t"
                 "n_present\tn_total\tbenefit\n")
        for s in synergies:
            slots = [[m] if isinstance(m, str) else list(m)
                     for m in s["requires"]]
            present_per_slot = [
                [m for m in slot if calls.get(m) in PRESENT_STATUSES]
                for slot in slots
            ]
            forbids = s.get("forbids", []) or []
            forbids_violated = [m for m in forbids
                                if calls.get(m) in PRESENT_STATUSES]
            forbids_satisfied = [m for m in forbids if m not in forbids_violated]

            n_slots_filled = sum(1 for p in present_per_slot if p)
            n_present = n_slots_filled + len(forbids_satisfied)
            n_total = len(slots) + len(forbids)
            comp = n_present / n_total if n_total else 0.0

            present_members = [m for p in present_per_slot for m in p]
            missing_slots = ["|".join(slot)
                             for slot, p in zip(slots, present_per_slot) if not p]

            fh.write("\t".join([
                s["name"], f"{comp:.3f}",
                status_for(comp, len(forbids_violated),
                           n_slots_filled, len(slots),
                           bool(s.get("single_organism"))),
                ",".join(present_members), ",".join(missing_slots),
                ",".join(forbids_satisfied), ",".join(forbids_violated),
                str(n_present), str(n_total),
                s.get("benefit", ""),
            ]) + "\n")


if __name__ == "__main__":
    main()
