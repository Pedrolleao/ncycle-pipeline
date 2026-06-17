#!/usr/bin/env python3
"""
build_custom_hmms.py — seed-curated HMM builder for Tier 2 / Tier 3 custom
HMM targets (custom_hmm: true in config/targets.yaml).

Single `build --target {id}` sub-command:

  build --target {id}
    Read targets/{id}/manifest.yaml -> expanded_candidates rows where keep:true
    Auto-fetch any missing accessions from UniProt (so the manifest is the
    single source of truth — curators can add accessions directly without
    re-running an expand step)
    Materialize refs.fasta (target_id||acc retag)
    CD-HIT cluster at cdhit_threshold (manifest override; default 0.7)
    MAFFT --auto -> refs.aln
    hmmbuild --name {target_id} -> {target_id}.hmm
    Update manifest.yaml > n_seqs_final + curation_date

The legacy phmmer-vs-UniRef90 `expand` sub-command was removed in P5.4.3
(2026-05-30) as dead code: recent batches (amoA_gamma, nosZ_clade2,
P5.2.1 N. multiformis/communis, P5.2.4 Chloroflexota) all curated their
expanded_candidates manually with keep:true and never invoked expand;
the UniRef90 download (~25 GB compressed / ~75 GB uncompressed) was
never triggered. cmd_build's auto-fetch of missing accessions makes
expand structurally unnecessary.

Reuses the UniProt fetch + header-retag pattern from build_blast_db.py.

Usage:
  python workflow/scripts/build_custom_hmms.py build --target nosZ_clade2
"""

from __future__ import annotations

import argparse
import datetime as _dt
import subprocess
import sys
from pathlib import Path

import yaml

from _common import ROOT
from build_blast_db import fetch_uniprot_fasta, retag_headers

TARGETS_DIR = ROOT / "targets"

MIN_SEQS = 8   # floor from WORK_PLAN section 4


def _target_dir(target_id: str) -> Path:
    d = TARGETS_DIR / target_id
    if not d.exists():
        sys.exit(f"error: {d} does not exist — create it and seed manifest.yaml first")
    return d


def _load_manifest(target_dir: Path) -> dict:
    mf = target_dir / "manifest.yaml"
    if not mf.exists():
        sys.exit(f"error: {mf} missing")
    with open(mf) as fh:
        return yaml.safe_load(fh) or {}


def _save_manifest(target_dir: Path, data: dict) -> None:
    mf = target_dir / "manifest.yaml"
    with open(mf, "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False)


# ── parsing helpers ─────────────────────────────────────────────────────────

def _bare_acc(tname: str) -> str:
    """Normalise a FASTA target name to a bare accession.

    Accepts SwissProt (`sp|Q9X4K0|XX_YYY`), TrEMBL (`tr|...|...`),
    UniRef (`UniRef90_A0A123`), retagged (`tid||acc`), or plain accessions.
    """
    if "||" in tname:
        return tname.split("||", 1)[1].split()[0]
    if "|" in tname:
        parts = tname.split("|")
        return parts[1] if len(parts) >= 2 else parts[0]
    if tname.startswith(("UniRef90_", "UniRef100_", "UniRef50_")):
        return tname.split("_", 1)[1]
    return tname.split()[0]


def _parse_fasta_str(text: str) -> dict[str, tuple[str, str]]:
    """Return {bare_acc: (header_line, sequence)} for a FASTA string."""
    out: dict[str, tuple[str, str]] = {}
    acc = None
    header = ""
    seq_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                out[acc] = (header, "".join(seq_lines))
            header = line.rstrip("\n")
            acc = _bare_acc(header[1:].split(None, 1)[0])
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if acc is not None:
        out[acc] = (header, "".join(seq_lines))
    return out


def _parse_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Return {bare_acc: (header_line, sequence)} for a FASTA file."""
    with open(path) as fh:
        return _parse_fasta_str(fh.read())


# ── build sub-command ───────────────────────────────────────────────────────

def cmd_build(args: argparse.Namespace) -> None:
    target = args.target
    tdir = _target_dir(target)
    manifest = _load_manifest(tdir)
    cands = manifest.get("expanded_candidates") or []
    kept = [c for c in cands if c.get("keep") is True]
    # MIN_SEQS default 8, override via manifest.min_seqs_override when the
    # family is biologically narrow (rare positives + extensive negatives is
    # still a working HMM). The override must be documented in manifest.notes.
    min_floor = int(manifest.get("min_seqs_override") or MIN_SEQS)
    if len(kept) < min_floor:
        sys.exit(f"error: only {len(kept)} sequence(s) marked keep: true "
                 f"(floor is {min_floor}). Either curate more rows or move "
                 f"this target to Tier 3 (BLAST-only).")
    print(f"[build] {target}: {len(kept)} curated sequence(s) → CD-HIT → MAFFT → hmmbuild")

    # Materialize refs.fasta from expanded.fasta, subset to kept accs.
    # Curators occasionally add accessions they found via fresh UniProt query
    # rather than from the expand step's output (happens when the seed was
    # wrong and the expansion was poisoned — e.g. hcnA's Q9I3F9 being an
    # alpha/beta hydrolase, not hcnA). Auto-fetch missing accessions from
    # UniProt so the manifest is the source of truth, not expanded.fasta.
    exp_records = _parse_fasta(tdir / "expanded.fasta")
    missing_accs = [c["acc"] for c in kept if c["acc"] not in exp_records]
    if missing_accs:
        print(f"[build] {target}: {len(missing_accs)} kept accession(s) not in "
              f"expanded.fasta — fetching from UniProt", file=sys.stderr)
        fresh = fetch_uniprot_fasta(missing_accs)
        if fresh.strip():
            fresh_records = _parse_fasta_str(fresh)
            exp_records.update(fresh_records)
    refs_chunks = []
    for c in kept:
        rec = exp_records.get(c["acc"])
        if rec is None:
            print(f"  ! kept accession {c['acc']} not in expanded.fasta and "
                  f"not fetchable from UniProt — skipping", file=sys.stderr)
            continue
        refs_chunks.append(f">{target}||{c['acc']}\n{rec[1]}\n")
    refs_fa = tdir / "refs.fasta"
    refs_fa.write_text("".join(refs_chunks))

    # CD-HIT redundancy clustering. Default 0.7 per WORK_PLAN §4, but highly
    # conserved narrow families (e.g. rusticyanin) collapse below the MIN_SEQS
    # floor at 0.7 — raise per-target via manifest.cdhit_threshold.
    cdhit_c = float(manifest.get("cdhit_threshold") or 0.7)
    cdhit_n = 5 if cdhit_c >= 0.7 else (4 if cdhit_c >= 0.6 else 3)
    cdhit_out = tdir / "refs.cdhit.fasta"
    subprocess.run(
        ["cd-hit", "-i", str(refs_fa), "-o", str(cdhit_out),
         "-c", f"{cdhit_c}", "-n", str(cdhit_n),
         "-T", str(args.threads), "-M", "0"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    # Replace refs.fasta with the clustered version.
    refs_fa.write_text(cdhit_out.read_text())
    cdhit_out.unlink(missing_ok=True)
    (tdir / "refs.fasta.clstr").unlink(missing_ok=True)
    Path(str(cdhit_out) + ".clstr").unlink(missing_ok=True)
    n_after = sum(1 for line in refs_fa.read_text().splitlines() if line.startswith(">"))
    if n_after < min_floor:
        sys.exit(f"error: only {n_after} sequence(s) after CD-HIT clustering "
                 f"(floor is {min_floor}). Curate more diverse sequences "
                 f"or raise cdhit_threshold toward 0.95.")

    # MAFFT alignment. --anysymbol handles selenocysteine (U) and pyrrolysine
    # (O), which appear in legitimate sequences like fdhF (E. coli formate
    # dehydrogenase F has Sec at residue 140).
    aln = tdir / "refs.aln"
    print(f"[build] {target}: aligning {n_after} sequence(s) with MAFFT")
    with open(aln, "w") as out:
        subprocess.run(
            ["mafft", "--anysymbol", "--auto", "--thread", str(args.threads), str(refs_fa)],
            check=True,
            stdout=out,
            stderr=subprocess.DEVNULL,
        )

    # hmmbuild.
    hmm = tdir / f"{target}.hmm"
    subprocess.run(
        ["hmmbuild", "--cpu", str(args.threads), "-n", target,
         str(hmm), str(aln)],
        check=True,
        stdout=subprocess.DEVNULL,
    )

    manifest["n_seqs_final"] = n_after
    manifest["curation_date"] = _dt.date.today().isoformat()
    _save_manifest(tdir, manifest)
    print(f"[build] {target}: wrote {hmm} ({n_after} sequence(s))")


# ── main ────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build", help="CD-HIT + MAFFT + hmmbuild from kept rows")
    p_build.add_argument("--target", required=True)
    p_build.add_argument("--threads", type=int, default=4)
    p_build.set_defaults(func=cmd_build)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
