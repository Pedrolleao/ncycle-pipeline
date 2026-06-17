#!/usr/bin/env python3
"""
detect_seed_leakage.py — flag hold-out (genome, target) cells whose BLAST-gate
seed set contains the SAME organism (or a congener) as the hold-out genome.

Motivation (2026-06-10 stats audit): the hold-out is meant to measure
generalization to genera ABSENT from training. But several BLAST-gate seeds were
curated from organisms that are conspecific (or congeneric) with hold-out
genomes — e.g. amoA_gamma seed D5BWX5 IS *Nitrosococcus halophilus* Nc4, which is
hold-out `Nhalophilus_Nc4`. On such a cell the BLAST gate trivially passes (the
query protein is ~identical to a seed), so the cell measures memorization, not
generalization, and should be excluded from (or sensitivity-tested against) the
hold-out score.

Method (reproducible, not hardcoded):
  1. Parse the per-target seed organisms from the `OS=` fields of EVERY active
     detector source: both BLAST-gate DIAMOND DBs (blast_gated_refs + unstable_refs;
     targets route to one or the other, config.yaml `databases:`) AND the custom-HMM
     training FASTAs (targets/<id>/expanded.fasta) — a seed leaks whether it sits in
     the BLAST gate or the HMM it trained.
  2. For each hold-out genome (taxonomy from HOLDOUT_TAXON, sourced from
     build_ground_truth.py's documented identities), check every ground-truth
     cell whose `expected == present` (only positives can be leakage-inflated):
       - SPECIES match (genus+species identical)  -> severity "species" (strong)
       - GENUS match  (genus identical)            -> severity "genus"   (moderate)
  3. Emit validation/holdout_seed_leakage.tsv (genome, target, severity, seed_acc,
     seed_organism). score_ncycle.py reads it to report a de-leaked hold-out.

Only `present` cells are flagged: a leaked seed can manufacture a true positive,
but cannot create a false positive on a GT-absent cell, so absent cells are clean.
"""
from __future__ import annotations
import csv, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Every active source from which a target's detector can "memorize" a seed:
# both routed BLAST-gate DBs + each custom-HMM training set.
SEED_SOURCES = [
    ROOT / "resources" / "blast_db" / "blast_gated_refs.fasta",
    ROOT / "resources" / "blast_db" / "unstable_refs.fasta",
    *sorted((ROOT / "targets").glob("*/expanded.fasta")),
]
GT = ROOT / "validation" / "ground_truth.tsv"
OUT = ROOT / "validation" / "holdout_seed_leakage.tsv"

# Hold-out genome -> (genus, species) — the documented identities from
# build_ground_truth.py HOLDOUTS (each line there names the organism). This is
# inherent panel metadata; the MATCHING below is what makes the detector general.
HOLDOUT_TAXON = {
    "Bdiazoefficiens_USDA110": ("bradyrhizobium", "diazoefficiens"),
    "Rpalustris_CGA009":       ("rhodopseudomonas", "palustris"),
    "Adehalogenans_2CP1":      ("anaeromyxobacter", "dehalogenans"),
    "Nhalophilus_Nc4":         ("nitrosococcus", "halophilus"),
    "Nhollandica_Lb":          ("nitrolancea", "hollandica"),
    "Sbrodae":                 ("scalindua", "brodae"),
    "Sstutzeri_F2a":           ("stutzerimonas", "stutzeri"),
    "Mcapsulatus_Bath":        ("methylococcus", "capsulatus"),
}
# Genus synonyms: Stutzerimonas stutzeri was long classified as Pseudomonas
# stutzeri; seeds may carry either genus name. Treat as the same organism.
GENUS_SYNONYMS = {"pseudomonas": "stutzerimonas"}  # only when species == stutzeri


def parse_seed_organisms() -> dict[str, list[tuple[str, str, str]]]:
    """target -> list of (genus, species, acc) parsed from the OS= fields of every
    active seed source (both gate DBs + custom-HMM training FASTAs). Headers in all
    sources share the `>target||acc ... OS=Genus species ... OX=` convention."""
    out: dict[str, list[tuple[str, str, str]]] = {}
    os_re = re.compile(r"OS=(.+?)\s+OX=")
    for src in SEED_SOURCES:
        if not src.exists():
            continue
        for line in src.read_text().splitlines():
            if not line.startswith(">") or "||" not in line:
                continue
            head = line[1:]
            target, _, rest = head.partition("||")
            acc = rest.split()[0] if rest else "?"
            m = os_re.search(head)
            if not m:
                continue
            toks = m.group(1).split()
            if len(toks) < 2:
                continue
            genus, species = toks[0].lower(), toks[1].lower().rstrip(".,")
            if species in ("sp", "sp.", "bacterium"):
                species = ""   # genus-only seed (e.g. "Bradyrhizobium sp.")
            out.setdefault(target, []).append((genus, species, acc))
    return out


def main() -> int:
    seeds = parse_seed_organisms()
    # hold-out positive cells
    holdout_pos: list[tuple[str, str]] = []
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("source") == "holdout_v3" and r["expected"] == "present":
                holdout_pos.append((r["genome"], r["target"]))

    rows = []
    for genome, target in holdout_pos:
        tax = HOLDOUT_TAXON.get(genome)
        if not tax:
            continue
        g_genus, g_species = tax
        best = None  # (severity_rank, severity, acc, organism)
        for (s_genus, s_species, acc) in seeds.get(target, []):
            s_genus = GENUS_SYNONYMS.get(s_genus, s_genus) if s_species == "stutzeri" else s_genus
            if s_genus != g_genus:
                continue
            if s_species and s_species == g_species:
                cand = (2, "species", acc, f"{s_genus} {s_species}")
            else:
                cand = (1, "genus", acc, f"{s_genus} {s_species or 'sp.'}")
            if best is None or cand[0] > best[0]:
                best = cand
        if best:
            rows.append({"genome": genome, "target": target, "severity": best[1],
                         "seed_acc": best[2], "seed_organism": best[3]})

    rows.sort(key=lambda r: (r["genome"], r["target"]))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "target", "severity",
                                           "seed_acc", "seed_organism"], delimiter="\t")
        w.writeheader(); w.writerows(rows)

    n_sp = sum(1 for r in rows if r["severity"] == "species")
    n_ge = sum(1 for r in rows if r["severity"] == "genus")
    print(f"[leakage] {len(rows)} leaked hold-out positive cells "
          f"({n_sp} species-level, {n_ge} genus-level) -> {OUT}")
    for r in rows:
        print(f"   {r['genome']:<24} {r['target']:<14} {r['severity']:<8} "
              f"{r['seed_acc']:<12} {r['seed_organism']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
