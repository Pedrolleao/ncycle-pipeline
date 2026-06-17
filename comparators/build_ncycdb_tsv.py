#!/usr/bin/env python3
"""
build_ncycdb_tsv.py — run NCycDB (Tu et al. 2019) on the training panel and export the
normalized TSV the benchmark harness consumes (columns: sample, family, count).

Faithful to NCyc's own profiler (NCycProfiler.PL), protein mode, DEFAULT settings:
    diamond makedb --in data/NCyc_100.faa --db data/NCyc_100      (done once)
    diamond blastp -k 1 -e 0.0001 -d data/NCyc_100 -q <proteome> -o <hits>
then best-hit subject id -> gene family via data/id2gene.map; a family's `count` for a
genome = number of query proteins whose top hit lands in that family (present iff > 0).

Training panel only (hold-outs excluded, matching compare_kofam.py / load_truth).
NCyc provenance: github.com/qichao1984/NCyc, NCyc_100.faa (219,146 seqs, 100% id), id2gene.map.
"""
from __future__ import annotations
import csv, subprocess, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NCYC = HERE / "NCyc" / "data"
DB = NCYC / "NCyc_100"
ID2GENE = NCYC / "id2gene.map"
ROOT = HERE.parent                       # ncycle-pipeline/
GT = ROOT / "validation" / "ground_truth.tsv"
PANEL = ROOT.parent / "test_panel"       # Nitrogen_Cycle/test_panel/*.faa
OUTDIR = HERE / "ncycdb_out"
OUT_TSV = ROOT / "validation" / "benchmark" / "ncycdb.tsv"
HOLDOUT_TAG = "holdout_v3"
EVALUE = "0.0001"                        # NCycProfiler.PL default
THREADS = "8"


def training_genomes() -> list[str]:
    train, hold = set(), set()
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            (hold if r.get("source") == HOLDOUT_TAG else train).add(r["genome"])
    return sorted(train - hold)


def load_id2gene() -> dict[str, str]:
    m = {}
    with open(ID2GENE) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                m[p[0]] = p[1]
    return m


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Run NCycDB on a proteome panel -> normalized TSV.")
    ap.add_argument("--proteome-dir", type=Path, default=PANEL,
                    help="dir of <genome>.faa (default: curated training panel)")
    ap.add_argument("--out", type=Path, default=OUT_TSV, help="output normalized TSV")
    ap.add_argument("--outdir", type=Path, default=OUTDIR, help="per-genome DIAMOND-hits cache dir")
    args = ap.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    if args.proteome_dir == PANEL:
        genomes = training_genomes()
    else:
        genomes = sorted(p.stem for p in args.proteome_dir.glob("*.faa"))
    id2gene = load_id2gene()
    print(f"[ncycdb] {len(genomes)} genomes from {args.proteome_dir}; {len(id2gene)} id->gene entries", file=sys.stderr)

    rows = []
    for g in genomes:
        faa = args.proteome_dir / f"{g}.faa"
        if not faa.exists():
            print(f"[ncycdb] WARN missing proteome {faa}", file=sys.stderr)
            continue
        hits = outdir / f"{g}.ncyc.tsv"
        if not hits.exists():
            cmd = ["diamond", "blastp", "-k", "1", "-e", EVALUE, "-p", THREADS,
                   "-d", str(DB), "-q", str(faa), "-o", str(hits), "--quiet"]
            subprocess.run(cmd, check=True)
        fam_count: dict[str, int] = defaultdict(int)
        seen_q: set[str] = set()
        with open(hits) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 2:
                    continue
                q, subj = c[0], c[1]
                if q in seen_q:        # -k 1 already keeps best hit; guard anyway
                    continue
                seen_q.add(q)
                fam = id2gene.get(subj)
                if fam:
                    fam_count[fam] += 1
        for fam, n in sorted(fam_count.items()):
            rows.append({"sample": g, "family": fam, "count": n})
        print(f"[ncycdb] {g}: {sum(fam_count.values())} mapped hits, {len(fam_count)} families", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "family", "count"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[ncycdb] wrote {args.out} ({len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
