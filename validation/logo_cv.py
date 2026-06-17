#!/usr/bin/env python3
"""
logo_cv.py — leave-one-clade-out cross-validation for a custom HMM (Frame C, generalization).

Ported from scycle-pipeline/validation/logo_cv.py to harmonize the two pipelines' reporting.

Question: does a custom HMM detect its target by capturing the conserved fold (generalizes to
clades NOT in its training refs), or is detection circular (it only works because the test
genome's homolog was a training ref)?

Method: read targets/{id}/manifest.yaml kept refs + their organism (ncycle stores organism text
in the `note` field; `organism` is also honoured if present), drop every ref whose organism text
matches a held-out clade keyword, rebuild the HMM (MAFFT --auto -> hmmbuild) from the remaining
refs, then hmmsearch the rebuilt HMM against each test proteome and report the best full-sequence
bitscore vs the calibrated TC. A held-out positive that still scores well above TC (and above the
negatives) is genuine fold-based generalization.

The N-cycle custom HMMs worth LOGO-testing are the homology-trap discriminators:
  nxrA / nxrB  (NOB clade vs denitrifier NarG/NarH)
  amoA         (AOB + comammox clade vs methanotroph PmoA)
  amoA_gamma   (gamma-AOB clade)
  nosZ_clade2  (clade-II atypical NosZ)

Usage examples:
  # Does the NOB nxrA HMM still fire on Nitrobacter when all Nitrobacter refs are held out?
  python validation/logo_cv.py --target nxrA \
      --holdout Nitrobacter \
      --test Nwinogradskyi_Nb255 Ngracilis_3211 Ninopinata_comammox \
             Paeruginosa_PAO1 Ecoli_K12_MG1655 Bsubtilis_168
  # Does the AOB amoA HMM generalize to comammox Nitrospira with all Nitrospira refs held out?
  python validation/logo_cv.py --target amoA --holdout Nitrospira \
      --test Ninopinata_comammox Nnitrosa_comammox Neuropaea_ATCC19718 Mfumariolicum_SolV
  (proteomes resolved from --panel, default ../test_panel; needs mafft + hmmbuild + hmmsearch)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def kept_refs(target: str) -> dict[str, str]:
    """acc -> organism text. ncycle keeps organism in `note`; honour `organism` too."""
    mf = yaml.safe_load(open(ROOT / "targets" / target / "manifest.yaml"))
    out = {}
    for c in mf.get("expanded_candidates", []):
        if c.get("keep") is True:
            out[c["acc"]] = f"{c.get('note', '')} {c.get('organism', '')}".strip()
    return out


def tc_of(target: str) -> float:
    """tc_cutoffs.tsv columns: profile_id<TAB>ga<TAB>tc. Custom HMMs use the TC column."""
    for line in (ROOT / "resources" / "hmm" / "tc_cutoffs.tsv").read_text().splitlines():
        p = line.split("\t")
        if p and p[0] == target and len(p) >= 3 and p[2]:
            return float(p[2])
    return float("nan")


def best_score(hmm: Path, faa: Path, work: Path) -> float | None:
    tbl = work / "out.tbl"
    subprocess.run(["hmmsearch", "--noali", "--tblout", str(tbl), str(hmm), str(faa)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    if not tbl.exists():
        return None
    for line in tbl.read_text().splitlines():
        if not line.startswith("#"):
            return float(line.split()[5])   # full-seq score column
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--holdout", nargs="+", required=True,
                    help="organism-name keywords; refs whose organism text matches any are dropped")
    ap.add_argument("--test", nargs="+", required=True, help="panel proteome stems to score")
    ap.add_argument("--panel", type=Path, default=ROOT.parent / "test_panel")
    args = ap.parse_args()

    refs = kept_refs(args.target)
    held = {a for a, o in refs.items() if any(k in o for k in args.holdout)}
    print(f"[logo] {args.target}: dropping {len(held)}/{len(refs)} refs matching {args.holdout}",
          file=sys.stderr)

    # subset refs.fasta to the non-held-out refs (header form: >target||ACC ...)
    seqs, acc = {}, None
    for line in (ROOT / "targets" / args.target / "refs.fasta").read_text().splitlines(keepends=True):
        if line.startswith(">"):
            acc = line[1:].split("||")[-1].split()[0].strip()
            seqs[acc] = [line]
        elif acc:
            seqs[acc].append(line)
    kept = [a for a in seqs if a not in held]

    tc = tc_of(args.target)
    with tempfile.TemporaryDirectory() as wd:
        work = Path(wd)
        (work / "na.fasta").write_text("".join(l for a in kept for l in seqs[a]))
        subprocess.run(f"mafft --auto {work}/na.fasta > {work}/na.aln",
                       shell=True, stderr=subprocess.DEVNULL, check=True)
        subprocess.run(["hmmbuild", "--amino", "-n", f"{args.target}_logo",
                        str(work / "logo.hmm"), str(work / "na.aln")],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print(f"\nleave-out {args.target} HMM ({len(kept)} refs) vs panel — TC={tc}")
        print(f"  {'proteome':<28}{'bitscore':>10}  call")
        for stem in args.test:
            faa = args.panel / f"{stem}.faa"
            s = best_score(work / "logo.hmm", faa, work) if faa.exists() else None
            call = "—" if s is None else (">=TC ✓" if s >= tc else "below TC")
            print(f"  {stem:<28}{(f'{s:.1f}' if s else 'none'):>10}  {call}")


if __name__ == "__main__":
    main()
