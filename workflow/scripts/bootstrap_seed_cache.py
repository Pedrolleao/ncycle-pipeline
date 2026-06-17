#!/usr/bin/env python3
"""
bootstrap_seed_cache.py — one-shot: pin the currently-built BLAST DB seeds
+ the currently-fetched Pfam HMM fallback into the cache (ROADMAP P5.1.5,
2026-05-30).

Walks the existing `resources/blast_db/{blast_gated_refs,unstable_refs}.fasta`
(both built from UniProt via build_blast_db.py), splits per accession, writes
each as `resources/.cache/seeds/<acc>.fasta` (canonical form: raw UniProt
header, no target retag), and records the manifest row with sha256 + length.

Also fingerprints the cached Pfam fallback HMMs at
`resources/.cache/ko_hmms/<pfam>.hmm` into `resources/.cache/pfam_release.txt`.

Designed to run once at the time of pinning. After this point, builds will
prefer the cache (P5.1.5 phase 2 — wired into build_blast_db.py and
build_hmm_db.py).

Usage:
  python workflow/scripts/bootstrap_seed_cache.py [--force]

  --force overwrites an existing manifest (audit the diff if you do this).
"""
from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

# Make sibling modules importable when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import ROOT, load_targets, blast_fallback_targets, blast_gated_targets  # noqa: E402
from _seed_cache import (  # noqa: E402
    SEED_CACHE_DIR, SEED_MANIFEST, PFAM_PIN,
    parse_fasta_by_accession, cache_seed, sha256_of, canonical_record,
    write_manifest, load_manifest,
)

BLAST_DIR = ROOT / "resources" / "blast_db"
PFAM_HMM_DIR = ROOT / "resources" / ".cache" / "ko_hmms"

GATED_FASTA = BLAST_DIR / "blast_gated_refs.fasta"
UNSTABLE_FASTA = BLAST_DIR / "unstable_refs.fasta"


def acc_to_target_from_retag(record: str) -> tuple[str | None, str | None]:
    """Read a retag header `>target||ACC sp|ACC|…` and return (target, ACC)."""
    line = record.splitlines()[0]
    if "||" not in line:
        return None, None
    tag = line[1:].split(" ", 1)[0]
    parts = tag.split("||", 1)
    if len(parts) != 2:
        return None, None
    return parts[0], parts[1]


def snapshot_blast_fasta(path: Path, fetch_date: str) -> list[dict]:
    """Read a build_blast_db.py output FASTA → cache rows."""
    if not path.exists():
        print(f"  ! {path} not present; skipping", file=sys.stderr)
        return []
    text = path.read_text()
    rows: list[dict] = []
    # parse_fasta_by_accession handles the retag form natively.
    for acc, record in parse_fasta_by_accession(text).items():
        target_id, acc2 = acc_to_target_from_retag(record)
        # Prefer the parsed accession; if retag agreed, both match.
        if not target_id:
            print(f"  ! no target_id retag for {acc} in {path.name}; "
                  f"skipping (cache only stores target-tagged seeds)",
                  file=sys.stderr)
            continue
        row = cache_seed(target_id, acc, record,
                         source=f"snapshot:{path.name}",
                         fetch_date=fetch_date)
        rows.append(row)
    return rows


def snapshot_pfams(fetch_date: str) -> list[dict]:
    """Hash each fetched Pfam HMM and write the pin file. Returns the rows
    for inclusion in pfam_release.txt."""
    pfams: list[dict] = []
    if not PFAM_HMM_DIR.exists():
        return pfams
    for hmm in sorted(PFAM_HMM_DIR.glob("PF*.hmm")):
        text = hmm.read_text()
        sha = sha256_of(canonical_record(text))
        pfams.append({
            "pfam_id": hmm.stem,
            "sha256": sha,
            "length_bytes": str(len(text)),
            "fetch_date": fetch_date,
            "source": "interpro:rest",
        })
    return pfams


def write_pfam_pin(rows: list[dict]) -> None:
    PFAM_PIN.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# resources/.cache/pfam_release.txt — Pfam fallback HMM pin (P5.1.5)",
             "# Each row is one Pfam-A HMM cached by build_hmm_db.py from InterPro.",
             "# Columns: pfam_id\\tsha256\\tlength_bytes\\tfetch_date\\tsource",
             ""]
    for r in rows:
        lines.append("\t".join([r["pfam_id"], r["sha256"], r["length_bytes"],
                                r["fetch_date"], r["source"]]))
    PFAM_PIN.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--force", action="store_true",
                    help="Overwrite an existing seeds_manifest.tsv")
    args = ap.parse_args()

    if SEED_MANIFEST.exists() and not args.force:
        existing = load_manifest()
        print(f"[bootstrap_seed_cache] manifest already exists with "
              f"{len(existing)} entries at {SEED_MANIFEST}")
        print("[bootstrap_seed_cache] pass --force to overwrite "
              "(audit the diff first)")
        return 0

    today = date.today().isoformat()
    by_key: dict[tuple[str, str], dict] = {}
    for fa in (GATED_FASTA, UNSTABLE_FASTA):
        print(f"[bootstrap_seed_cache] snapshotting {fa.relative_to(ROOT)}…")
        new_rows = snapshot_blast_fasta(fa, today)
        print(f"  + cached {len(new_rows)} accession(s)")
        for r in new_rows:
            by_key[(r["target_id"], r["accession"])] = r
    rows = list(by_key.values())

    write_manifest(rows)
    print(f"[bootstrap_seed_cache] wrote {SEED_MANIFEST} "
          f"({len(rows)} accessions across "
          f"{len({r['target_id'] for r in rows})} targets)")
    print(f"[bootstrap_seed_cache] per-accession snapshots in {SEED_CACHE_DIR}")

    pfam_rows = snapshot_pfams(today)
    if pfam_rows:
        write_pfam_pin(pfam_rows)
        print(f"[bootstrap_seed_cache] wrote {PFAM_PIN} "
              f"({len(pfam_rows)} Pfam HMM(s) pinned)")
    else:
        print(f"[bootstrap_seed_cache] no Pfam HMMs found under "
              f"{PFAM_HMM_DIR.relative_to(ROOT)}; skipping Pfam pin")

    # Sanity: every target with blast_refs_uniprot should have at least one cached seed.
    targets = load_targets()["targets"]
    cached_by_target: dict[str, list[str]] = {}
    for r in rows:
        cached_by_target.setdefault(r["target_id"], []).append(r["accession"])
    expected_targets = set()
    for t in [*blast_fallback_targets(targets), *blast_gated_targets(targets)]:
        if t.get("blast_refs_uniprot"):
            expected_targets.add(t["id"])
    missing = expected_targets - cached_by_target.keys()
    if missing:
        print(f"  ! {len(missing)} target(s) with blast_refs_uniprot but no "
              f"cached seeds: {sorted(missing)}", file=sys.stderr)
        return 1
    print(f"[bootstrap_seed_cache] verified: every target with "
          f"blast_refs_uniprot is represented in the cache")
    return 0


if __name__ == "__main__":
    sys.exit(main())
