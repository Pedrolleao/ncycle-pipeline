#!/usr/bin/env python3
"""
verify_seeds.py — read-only drift check for the UniProt seed cache + Pfam pin.

ROADMAP P5.1.5 (2026-05-30). The seed cache at `resources/.cache/seeds/`
+ pin at `resources/.cache/seeds_manifest.tsv` make BLAST DB rebuilds
reproducible (see `workflow/scripts/_seed_cache.py`). This script
re-fetches every cached seed accession from UniProt + every pinned Pfam
HMM from InterPro and diffs the sha256s against the pin. Output:

  - match     — UniProt content sha256 == pinned sha256  (no drift)
  - mismatch  — UniProt has revised the entry since pinning
  - missing   — UniProt no longer returns the accession (deprecated/merged)
  - error     — network / fetch failure

Exit code 0 if everything matches; 1 if any mismatch or missing; 2 on
network errors only.

Intended uses:
  - `make verify-seeds` periodic CI / pre-publication check
  - debugging: did UniProt update a seed underneath us?

This is read-only — does NOT update the cache or manifest. To intentionally
refresh after auditing the report, use:
  python workflow/scripts/build_blast_db.py --force --refresh-seeds
"""
from __future__ import annotations

import argparse
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from _seed_cache import (  # noqa: E402
    SEED_CACHE_DIR, SEED_MANIFEST, PFAM_PIN,
    parse_fasta_by_accession, canonical_record, sha256_of, load_manifest,
)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
INTERPRO_HMM_URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pf}/?annotation=hmm"


def fetch_uniprot_batch(accessions: list[str]) -> dict[str, str]:
    """Fetch FASTA for a batch of accessions and parse to {acc -> record}."""
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{a}" for a in accessions)
    params = urllib.parse.urlencode({"query": query, "format": "fasta",
                                      "size": "500"})
    url = f"{UNIPROT_BASE}/search?{params}"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode()
    except Exception as e:
        print(f"  ! UniProt fetch failed for batch of {len(accessions)}: {e}",
              file=sys.stderr)
        return {}
    return parse_fasta_by_accession(text)


def fetch_pfam_hmm(pf: str) -> str | None:
    """Return raw Pfam HMM text from InterPro, or None on failure."""
    import gzip
    url = INTERPRO_HMM_URL.format(pf=pf)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = r.read()
        try:
            return gzip.decompress(data).decode()
        except OSError:
            return data.decode(errors="replace")
    except Exception as e:
        print(f"  ! InterPro fetch failed for {pf}: {e}", file=sys.stderr)
        return None


def verify_uniprot(batch_size: int = 50) -> dict:
    """Compare cached UniProt seeds against live UniProt. Returns counts dict."""
    manifest = load_manifest()
    if not manifest:
        print("[verify_seeds] no manifest at "
              f"{SEED_MANIFEST.relative_to(ROOT)} — run "
              "workflow/scripts/bootstrap_seed_cache.py first",
              file=sys.stderr)
        return {"status": "no_manifest"}

    by_target = defaultdict(list)
    for r in manifest:
        by_target[r["target_id"]].append(r)

    n_match = n_mismatch = n_missing = n_error = 0
    mismatches: list[tuple[str, str]] = []
    missings:   list[tuple[str, str]] = []
    print(f"[verify_seeds] verifying {len(manifest)} pinned UniProt seeds "
          f"across {len(by_target)} target(s)…")
    # Fetch in batches of ~batch_size per target group to balance request count
    # vs URL length.
    for tid, rows in sorted(by_target.items()):
        accs = [r["accession"] for r in rows]
        all_records: dict[str, str] = {}
        for i in range(0, len(accs), batch_size):
            chunk = accs[i : i + batch_size]
            all_records.update(fetch_uniprot_batch(chunk))
            time.sleep(0.4)
        for r in rows:
            acc = r["accession"]
            pinned_sha = r["sha256"]
            live_record = all_records.get(acc)
            if live_record is None:
                n_missing += 1
                missings.append((tid, acc))
                continue
            live_sha = sha256_of(canonical_record(live_record))
            if live_sha == pinned_sha:
                n_match += 1
            else:
                n_mismatch += 1
                mismatches.append((tid, acc))
        print(f"  {tid:<16} {sum(1 for r in rows if all_records.get(r['accession']) and sha256_of(canonical_record(all_records[r['accession']])) == r['sha256']):>3}/{len(rows)} match")

    print()
    print(f"[verify_seeds] UniProt summary: {n_match} match, {n_mismatch} "
          f"mismatch (silent UniProt revision), {n_missing} missing "
          f"(deprecated/merged), {n_error} error")
    if mismatches:
        print("  mismatches (re-pin after auditing):")
        for tid, acc in mismatches:
            print(f"    {tid:<16} {acc}")
    if missings:
        print("  missings (UniProt no longer returns these):")
        for tid, acc in missings:
            print(f"    {tid:<16} {acc}")
    return {"status": "ok",
            "match": n_match, "mismatch": n_mismatch,
            "missing": n_missing, "error": n_error}


def verify_pfam() -> dict:
    """Compare cached Pfam HMMs against InterPro live. Returns counts dict."""
    if not PFAM_PIN.exists():
        print(f"[verify_seeds] no Pfam pin at {PFAM_PIN.relative_to(ROOT)}; "
              "skipping Pfam verification")
        return {"status": "no_pin"}
    pinned: dict[str, str] = {}
    for line in PFAM_PIN.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            pinned[parts[0]] = parts[1]
    if not pinned:
        return {"status": "empty_pin"}

    n_match = n_mismatch = n_error = 0
    mismatches: list[str] = []
    print(f"\n[verify_seeds] verifying {len(pinned)} pinned Pfam HMM(s)…")
    for pf, pinned_sha in pinned.items():
        text = fetch_pfam_hmm(pf)
        if text is None:
            n_error += 1
            print(f"  {pf} ERROR (fetch failed)")
            continue
        live_sha = sha256_of(canonical_record(text))
        if live_sha == pinned_sha:
            n_match += 1
            print(f"  {pf} match")
        else:
            n_mismatch += 1
            mismatches.append(pf)
            print(f"  {pf} MISMATCH (pin={pinned_sha[:12]}…, live={live_sha[:12]}…)")
    print(f"[verify_seeds] Pfam summary: {n_match} match, {n_mismatch} "
          f"mismatch, {n_error} error")
    return {"status": "ok", "match": n_match,
            "mismatch": n_mismatch, "error": n_error}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--skip-uniprot", action="store_true",
                    help="Skip UniProt verification (Pfam only)")
    ap.add_argument("--skip-pfam", action="store_true",
                    help="Skip Pfam verification (UniProt only)")
    args = ap.parse_args()

    uni = pfam = None
    if not args.skip_uniprot:
        uni = verify_uniprot()
    if not args.skip_pfam:
        pfam = verify_pfam()

    drift = 0
    error = 0
    for d in (uni, pfam):
        if d is None or d.get("status") not in ("ok",):
            continue
        drift += d.get("mismatch", 0) + d.get("missing", 0)
        error += d.get("error", 0)
    if drift > 0:
        print(f"\n[verify_seeds] DRIFT DETECTED: {drift} seed(s) have changed "
              "upstream since pinning. Re-pin via "
              "`python workflow/scripts/build_blast_db.py --force --refresh-seeds` "
              "after auditing what changed.")
        return 1
    if error > 0:
        print(f"\n[verify_seeds] {error} fetch error(s) — network issue; "
              "re-run later.")
        return 2
    print("\n[verify_seeds] OK: all pinned seeds match upstream.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
