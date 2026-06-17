#!/usr/bin/env python3
"""
build_metabolic_tsv.py — normalize METABOLIC v4.0 (Zhou et al. 2022) output → metabolic.tsv.

Output columns: genome, ko   (one present KO per row). A KO is present in a genome iff
METABOLIC's per-genome KEGG result lists a protein hit for it (KOfam adaptive thresholds,
applied by METABOLIC). Mapped to ncycle targets via the same KO->target table as the raw
KofamScan baseline (adapters.load_metabolic) — so METABOLIC's trap behavior is governed by the
same shared-KO ambiguity as KofamScan, by construction.

Source: METABOLIC_out/KEGG_identifier_result/<genome>.result.txt — TSV `KO <tab> count`
(col2 non-empty => present). Run AFTER METABOLIC-G.pl completes.
Note: METABOLIC's bundled KOfam version lacks K10534 (NR) and K17877 (nasD) — those two
ncycle targets are KO-version coverage gaps for METABOLIC (documented, <=2 cells).
"""
from __future__ import annotations
import csv, glob, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KEGG_DIR = HERE / "metabolic_out" / "KEGG_identifier_result"
OUT_TSV = HERE.parent / "validation" / "benchmark" / "metabolic.tsv"


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Normalize METABOLIC KEGG_identifier_result -> (genome, ko) TSV.")
    ap.add_argument("--kegg-dir", type=Path, default=KEGG_DIR,
                    help="METABOLIC_out/KEGG_identifier_result dir")
    ap.add_argument("--out", type=Path, default=OUT_TSV, help="output normalized TSV")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.kegg_dir / "*.result.txt")))
    if not files:
        print(f"[metabolic] no *.result.txt in {args.kegg_dir}", file=sys.stderr)
        return 1
    rows = []
    for fp in files:
        g = Path(fp).name[: -len(".result.txt")]
        n = 0
        with open(fp) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0].startswith("K") and p[1].strip():
                    rows.append({"genome": g, "ko": p[0]})
                    n += 1
        print(f"[metabolic] {g}: {n} present KOs", file=sys.stderr)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "ko"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"[metabolic] wrote {args.out} ({len(rows)} rows, {len(files)} genomes)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
