#!/usr/bin/env python3
"""
filter_proteome_by_hmm.py — produce a HMM-hit-filtered proteome FASTA for
DIAMOND-blastp (ROADMAP P5.4.2, optional efficiency path, 2026-05-31).

The DIAMOND-blastp rules in `workflow/rules/protein_mode.smk` currently
query the full proteome FASTA against the small target seed database. On
typical proteomes (3-5k proteins), >95% of proteins have zero HMM hits
across all KO + Pfam + custom profiles — DIAMOND aligns those to no
purpose, wasting most of the BLAST step's runtime. For 1000-genome
studies, filtering the proteome to only HMM-hit proteins before DIAMOND
cuts BLAST cost ~10× without changing any prediction.

This script is the optional pre-filter: read a proteome + hmmscan tblout,
emit a filtered FASTA containing only the proteins with at least one
HMM hit at the per-target threshold.

  ⚠️ NOT PREDICTION-NEUTRAL for ncycle (verified 2026-06-09 by wiring it in
  and re-scoring). The BLAST gate is COMPLEMENTARY to the HMMs, not nested
  under them: some targets are detected by curated BLAST seeds on divergent
  proteins that produce NO hmmscan hit at all (e.g. assimilatory narB/K00367
  in E. coli, nasD/K17877 in N. europaea — their proteins are absent from the
  hmmscan tblout entirely, yet BLAST matches them at high identity). Filtering
  the DIAMOND query by hmmscan hits drops those proteins → those calls flip to
  absent → training micro-F1 1.00 → 0.996 (2 FN), the regression gate FAILS.
  The wiring exists (rule `filter_proteome_by_hmm` + config
  options.filter_blast_query_by_hmm) but the toggle DEFAULTS OFF. Only enable
  it for throughput experiments where a small recall loss on BLAST-only targets
  is acceptable — never for an accuracy-critical run. (It also only speeds the
  already-cheap DIAMOND-vs-seeds step, not the hmmscan bottleneck that dominates
  large runs, so it is not the right lever for GTDB-scale efficiency.)

Usage (CLI):
  python workflow/scripts/filter_proteome_by_hmm.py \\
      --proteome results/SAMPLE/prodigal/SAMPLE.faa \\
      --hmmscan results/SAMPLE/hmm/SAMPLE.hmmscan.tsv \\
      --out results/SAMPLE/prodigal/SAMPLE.hmm_filtered.faa

Wired into the Snakemake DAG (2026-06-09) behind config
`options.filter_blast_query_by_hmm` (rule `filter_proteome_by_hmm` →
`diamond_blastp_unstable`/`_gated` query the filtered FASTA; hmmscan +
apply_rules gene-coordinate synteny always use the FULL proteome). The toggle
DEFAULTS OFF for the reason in the warning above (regression FAILS when on).

Filtering policy: a protein is kept if it has AT LEAST ONE hmmscan row,
regardless of bitscore. This is permissive on purpose — TC-filtering happens
later in apply_rules. NOTE the policy assumes BLAST-relevant proteins also have
an hmmscan row; that assumption does NOT hold for ncycle's BLAST-only targets
(see warning), which is why enabling the filter is lossy here.

Empirical retention rate (2026-05-31, 6 representative panel samples):

  Pdenitrificans_PD1222   78/5019  (1.6%)
  Adehalogenans_2CP1      88/4477  (2.0%)
  Sstutzeri_F2a           80/4195  (1.9%)
  Mcapsulatus_Bath        67/2971  (2.3%)
  Sbrodae                 73/3663  (2.0%)
  Nwinogradskyi_Nb255     41/3262  (1.3%)

Mean retention ~1.85% — DIAMOND would query ~80 proteins instead of ~4000,
a **~50× cost cut** (better than the audit's "~10×" estimate). For a
1000-genome run, this translates to ~2 orders of magnitude less DIAMOND
runtime without changing any prediction (the dropped proteins by
construction have zero HMM hits and therefore zero target-specific BLAST
evidence in apply_rules).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def proteins_with_hits(hmmscan_tsv: Path) -> set[str]:
    """Return the set of protein IDs that have at least one hmmscan hit.

    The pipeline writes `hmmscan --domtblout` per
    workflow/rules/protein_mode.smk; column 4 is the query name (the
    protein ID from the proteome FASTA), column 1 is the target (the HMM
    profile name).
    """
    hits: set[str] = set()
    with open(hmmscan_tsv) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            hits.add(parts[3])
    return hits


def filter_fasta(proteome: Path, keep_ids: set[str], out: Path) -> tuple[int, int]:
    """Stream-filter proteome.faa → out.faa keeping only records whose
    header ID is in keep_ids. Returns (n_total, n_kept)."""
    n_total = n_kept = 0
    keep = False
    with open(proteome) as src, open(out, "w") as dst:
        for line in src:
            if line.startswith(">"):
                n_total += 1
                # The protein ID is the first whitespace-separated token after `>`.
                pid = line[1:].split(None, 1)[0]
                keep = pid in keep_ids
                if keep:
                    n_kept += 1
                    dst.write(line)
            elif keep:
                dst.write(line)
    return n_total, n_kept


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--proteome", type=Path, required=True,
                    help="Input proteome FASTA (e.g. prodigal output)")
    ap.add_argument("--hmmscan", type=Path, required=True,
                    help="hmmscan --tblout TSV from the same proteome")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output filtered FASTA")
    args = ap.parse_args()

    if not args.proteome.is_file():
        sys.exit(f"error: proteome not found: {args.proteome}")
    if not args.hmmscan.is_file():
        sys.exit(f"error: hmmscan tsv not found: {args.hmmscan}")

    keep_ids = proteins_with_hits(args.hmmscan)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_total, n_kept = filter_fasta(args.proteome, keep_ids, args.out)
    pct = (100.0 * n_kept / n_total) if n_total else 0.0
    print(f"[filter_proteome_by_hmm] {args.proteome.name}: "
          f"kept {n_kept}/{n_total} ({pct:.1f}%) HMM-hit proteins → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
