#!/usr/bin/env python3
"""
build_dram_tsv.py — normalize DRAM (Shaffer et al. 2020) output → dram.tsv for the harness.

Output columns: genome, function, present   (function = DRAM KEGG-module label;
present = 1 iff DRAM detected >=1 of that module's member KOs in the genome — the
definition in adapters.load_dram). Faithful to the pre-registered DRAM_STEP_MAP, which
maps each module to the N-cycle step(s) it covers (coarse; expanded to subunits by the engine).

Inputs:
  dram_out/annotations.tsv                : per-gene calls; cols 'fasta' (genome) + 'ko_id'
  <DRAM_data>/forms/genome_summary_form*  : module -> member KO (gene_id col = KO, module col = label)
Run AFTER `DRAM.py annotate_genes` + `DRAM.py distill`.
"""
from __future__ import annotations
import csv, glob, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANNOT = HERE / "dram_out" / "annotations.tsv"
FORM_GLOB = "/home/dmin/Grants/Sulfur_Cycle/comparators/DRAM_data/forms/genome_summary_form*.tsv"
OUT_TSV = HERE.parent / "validation" / "benchmark" / "dram.tsv"

sys.path.insert(0, str(HERE.parent / "validation" / "benchmark"))
from adapters import DRAM_STEP_MAP  # the pre-registered module->step map (8 N modules)


def module_to_kos() -> dict[str, set[str]]:
    form = sorted(glob.glob(FORM_GLOB))[-1]
    want = {k.lower() for k in DRAM_STEP_MAP}
    m2k: dict[str, set[str]] = defaultdict(set)
    with open(form) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mod = (r.get("module") or "").strip()
            if mod.lower() in want and (r.get("gene_id") or "").startswith("K"):
                m2k[mod.lower()].add(r["gene_id"].strip())
    print(f"[dram] {form}: {len(m2k)} N-modules, "
          f"{sum(len(v) for v in m2k.values())} member KOs", file=sys.stderr)
    return m2k


def genome_kos(annot: Path = ANNOT) -> dict[str, set[str]]:
    g2k: dict[str, set[str]] = defaultdict(set)
    with open(annot) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        # genome name == 'fasta'; gene headers were prefixed '<genome>__' so strip that too
        for r in rdr:
            g = (r.get("fasta") or "").strip()
            ko = (r.get("ko_id") or "").strip()
            if g and ko and ko != "":
                g2k[g].add(ko)
    print(f"[dram] {len(g2k)} genomes in annotations.tsv", file=sys.stderr)
    return g2k


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Normalize DRAM annotations → dram.tsv")
    ap.add_argument("--annotations", type=Path, default=ANNOT)
    ap.add_argument("--out", type=Path, default=OUT_TSV)
    args = ap.parse_args()
    out_tsv = args.out
    m2k = module_to_kos()
    g2k = genome_kos(args.annotations)
    # canonical module label (original case) for each lowercased key, from DRAM_STEP_MAP order
    label = {}
    for k in DRAM_STEP_MAP:
        label[k.lower()] = k
    rows = []
    for g in sorted(g2k):
        for modlc, kos in m2k.items():
            present = 1 if (g2k[g] & kos) else 0
            rows.append({"genome": g, "function": label[modlc], "present": present})
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "function", "present"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    npos = sum(r["present"] for r in rows)
    print(f"[dram] wrote {out_tsv} ({len(rows)} rows, {npos} present)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
