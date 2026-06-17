"""
_seed_cache.py — shared helpers for the UniProt-seed snapshot cache
(ROADMAP P5.1.5, 2026-05-30).

The cache pins UniProt seed FASTAs to a local snapshot so BLAST DB rebuilds
are reproducible even if UniProt revises an accession (rare but happens
when sequences are updated, deprecated, or merged into another entry).
Layout:

  resources/.cache/seeds/<acc>.fasta     one file per UniProt accession
                                          (the raw UniProt FASTA record,
                                          no target_id retag)
  resources/.cache/seeds_manifest.tsv     columns: target_id, accession,
                                          sha256, length, fetch_date,
                                          source, description

A build script reads the cache first; if a needed accession is missing,
it falls back to a live UniProt fetch and writes the new snapshot.
`validation/verify_seeds.py` re-fetches everything and diffs against the
pinned sha256s to detect upstream revisions.

Helpers exposed:
  - SEED_CACHE_DIR / SEED_MANIFEST    path constants
  - parse_fasta_by_accession(text)     {acc -> single-record FASTA}
  - canonical_record(...)              normalize a single record for hashing
  - sha256_of(text)                    hex sha256
  - load_manifest()                    list[dict]
  - write_manifest(rows)               atomic write
  - cache_seed(target_id, acc, text)   write snapshot + return manifest row
  - read_cached_fasta(acc)             cached FASTA or None
"""
from __future__ import annotations

import hashlib
import re
from datetime import date
from pathlib import Path

from _common import ROOT

SEED_CACHE_DIR = ROOT / "resources" / ".cache" / "seeds"
SEED_MANIFEST  = ROOT / "resources" / ".cache" / "seeds_manifest.tsv"
PFAM_PIN       = ROOT / "resources" / ".cache" / "pfam_release.txt"
DEFAULT_SOURCE = "uniprot:rest"

# UniProt-style headers are >sp|ACC|NAME … or >tr|ACC|NAME … . Custom-prefixed
# headers (target_id||ACC …) also appear in some snapshot sources.
_ACC_RE = re.compile(r"^>(?:[^|]+\|\|)?(?:sp|tr)\|([A-Z0-9]+)\|", re.I)
_RETAG_ACC_RE = re.compile(r"^>(?:[^|]+\|\|)([A-Z0-9]+)\b", re.I)


def parse_fasta_by_accession(text: str) -> dict[str, str]:
    """Return {accession -> single-record FASTA (header + seq, newline-terminated)}.

    Handles both raw UniProt (>sp|ACC|…) and the build_blast_db.py retag style
    (>target||ACC sp|ACC|…). For retag style the accession before the space
    is preferred; for raw it's the UniProt accession in the second `|`-field.
    """
    out: dict[str, str] = {}
    current_acc: str | None = None
    current_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current_acc:
                out[current_acc] = "\n".join(current_lines) + "\n"
            current_lines = [line]
            m = _RETAG_ACC_RE.match(line) or _ACC_RE.match(line)
            current_acc = m.group(1) if m else None
        else:
            current_lines.append(line)
    if current_acc:
        out[current_acc] = "\n".join(current_lines) + "\n"
    return out


def strip_retag(record: str) -> str:
    """If a record's header is in retag form (>target||ACC sp|ACC|…),
    return it with the retag stripped so the canonical record is the raw
    UniProt FASTA. This makes the cache identical regardless of which
    build script wrote it."""
    lines = record.splitlines()
    if not lines:
        return record
    h = lines[0]
    if "||" in h and h.startswith(">"):
        # Drop "target_id||ACC " prefix; keep the rest of the header as-is.
        # Format example: ">amoA||Q04507 sp|Q04507|AMOA_NITEU Ammonia …"
        rest = h.split(" ", 1)
        if len(rest) == 2:
            tail = rest[1]
            if tail.startswith("sp|") or tail.startswith("tr|"):
                lines[0] = ">" + tail
    return "\n".join(lines) + "\n"


def canonical_record(record: str) -> str:
    """Normalize a single FASTA record for hashing: strip retag, normalize
    line endings, drop trailing blank lines."""
    record = strip_retag(record)
    text = "\n".join(line.rstrip() for line in record.splitlines() if line.strip())
    return text + "\n"


def sha256_of(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def cache_seed(target_id: str, accession: str, record: str,
               source: str = DEFAULT_SOURCE,
               fetch_date: str | None = None) -> dict:
    """Write a per-accession snapshot to the cache and return the manifest row.

    `record` is one FASTA record (header + seq, newline-terminated). It's
    normalized via `canonical_record()` before sha256+write so the on-disk
    file is the canonical form (raw UniProt header, no target retag).
    """
    SEED_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    canon = canonical_record(record)
    seq = "".join(l for l in canon.splitlines() if not l.startswith(">"))
    desc = canon.splitlines()[0] if canon else ""
    sha = sha256_of(canon)
    out = SEED_CACHE_DIR / f"{accession}.fasta"
    out.write_text(canon)
    return {
        "target_id": target_id,
        "accession": accession,
        "sha256": sha,
        "length": str(len(seq)),
        "fetch_date": fetch_date or date.today().isoformat(),
        "source": source,
        "description": desc[:120],
    }


def read_cached_fasta(accession: str) -> str | None:
    """Return cached FASTA for an accession, or None if not cached."""
    p = SEED_CACHE_DIR / f"{accession}.fasta"
    return p.read_text() if p.exists() else None


_MANIFEST_COLS = ["target_id", "accession", "sha256", "length",
                  "fetch_date", "source", "description"]


def load_manifest() -> list[dict]:
    if not SEED_MANIFEST.exists():
        return []
    out = []
    with open(SEED_MANIFEST) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            out.append(dict(zip(header, parts)))
    return out


def write_manifest(rows: list[dict]) -> None:
    SEED_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(SEED_MANIFEST, "w") as fh:
        fh.write("\t".join(_MANIFEST_COLS) + "\n")
        for r in sorted(rows, key=lambda x: (x.get("target_id", ""),
                                              x.get("accession", ""))):
            fh.write("\t".join(str(r.get(c, "")) for c in _MANIFEST_COLS) + "\n")
