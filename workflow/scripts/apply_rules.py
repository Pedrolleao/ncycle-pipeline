#!/usr/bin/env python3
"""
apply_rules.py — combine HMM (hmmsearch) + DIAMOND-blastp evidence, emit one row
per target into calls/ncycle_calls.tsv.

Detection model:
  - HMM evidence is an hmmsearch domtblout (KOfam KO profiles + custom clade HMMs
    + a Pfam fallback); see parse_hmmscan_domtbl.
  - A "BLAST hit" is a DIAMOND hit against a curated UniProt reference tagged for
    this target (header convention: >{target_id}||{acc}).
  - Status set: confirmed | domain-only | narrow-no-IPR | disqualified.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

KO_RE = re.compile(r"^K\d{4,}$")


# ───────────────────────────── parsers ───────────────────────────────────────

def parse_hmmscan_domtbl(path: Path, tc_cutoffs: dict[str, float],
                         qcov_cutoffs: dict[str, float] | None = None
                        ) -> dict[str, dict]:
    """Return {protein_id: {pfam: {pfam_id: best_evalue},
                             custom: {target_id: best_evalue}}}.

    Pfam profiles are identified by target_acc starting with "PF" (e.g.
    PF13435.10). Custom HMMs (built via build_custom_hmms.py with
    `hmmbuild --name target_id`) carry the target_id in target_name and
    typically have no ACC field (`-`).

    `qcov_cutoffs` is a per-target {target_id: min_query_coverage} dict.
    When set, domains with env-coverage of the query protein below the
    threshold are dropped — prevents small-domain HMMs (e.g. 104-aa hcnA)
    from false-positive-hitting the matching domain of a large multidomain
    protein (~950 aa Fe-S-containing oxidoreductase).
    """
    qcov_cutoffs = qcov_cutoffs or {}
    out: dict[str, dict] = defaultdict(lambda: {"pfam": {}, "custom": {}, "ko": {}})
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            # NOTE: produced by `hmmsearch --domtblout` (profiles vs the proteome),
            # which is several-fold faster than hmmscan for a small profile DB.
            # hmmsearch swaps the target/query roles vs hmmscan: the SEQUENCE
            # (protein) is the target (parts[0], tlen=parts[2]) and the PROFILE is
            # the query (parts[3], acc parts[4]). The env coords (parts[19/20]) are
            # on the sequence in both tools, and the full-seq score (parts[7]) is the
            # protein-vs-profile score in both — so only the name/acc/len columns move.
            prof_name = parts[3]                  # KO / custom-target / Pfam profile name
            target_acc = parts[4]                 # "PF#####.NN" or "-"
            qlen = int(parts[2])                  # protein (sequence) length, for qcov
            query_name = parts[0]                 # protein id
            target_name = prof_name
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            env_from = int(parts[19])
            env_to = int(parts[20])
            if target_acc.startswith("PF"):
                key = target_acc.split(".")[0]
                bucket = "pfam"
            elif KO_RE.match(target_name) or KO_RE.match(target_acc):
                # KOfam profile — name (or acc) is the KO number, e.g. K02588.
                key = target_name if KO_RE.match(target_name) else target_acc
                bucket = "ko"
            else:
                key = target_name
                bucket = "custom"
            tc = tc_cutoffs.get(key)
            # Apply TC cutoff if available, else trust the e-value filter that
            # hmmscan already applied via -E.
            if tc is not None and full_score < tc:
                continue
            # Apply per-target query-coverage filter (defends against
            # small-domain HMMs over-calling on multidomain proteins).
            min_qcov = qcov_cutoffs.get(key)
            if min_qcov is not None and qlen > 0:
                qcov = (env_to - env_from + 1) / qlen
                if qcov < min_qcov:
                    continue
            entry = out[query_name][bucket]
            prev = entry.get(key)
            if prev is None or full_evalue < prev:
                entry[key] = full_evalue
    return out


def parse_blast_tsv(path: Path) -> dict[str, dict[str, list[tuple[str, float, float]]]]:
    """Return {protein_id: {target_id: [(uniprot_acc, pident, evalue), …]}}.
    Subject header convention: target_id||uniprot_acc."""
    out: dict[str, dict] = defaultdict(lambda: defaultdict(list))
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qid = parts[0]
            sid = parts[1]
            if "||" not in sid:
                continue
            target_id, acc = sid.split("||", 1)
            pident = float(parts[2])
            evalue = float(parts[10])
            out[qid][target_id].append((acc.split(" ", 1)[0], pident, evalue))
    return out


def load_tc_cutoffs(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2]:
                try:
                    out[parts[0]] = float(parts[2])
                except ValueError:
                    pass
    return out


# ───────────────────────────── per-target evaluation ─────────────────────────

def evaluate_target(target: dict,
                    pfam_hits: dict[str, dict[str, float]],
                    custom_hits: dict[str, dict[str, float]],
                    ko_hits: dict[str, dict[str, float]],
                    blast_hits: dict[str, dict],
                    pident_default: float, qcov_default: float,
                    ) -> tuple[str, dict] | None:
    """Return (status, evidence_dict) for the BEST protein matching this target,
    or None if nothing matches.

    Status set:
      confirmed     — sequence signature met (Pfam OR custom HMM) AND BLAST hit
                       (or signature-only when no fallback is configured)
      domain-only   — signature met, no BLAST hit (only meaningful when fallback set)
      narrow-no-IPR — BLAST hit, signature not met
      disqualified  — requires_blast_for_confirmation set, signature met, BLAST missing
    """
    tid = target["id"]
    is_custom = bool(target.get("custom_hmm"))
    req_pfams = set(target.get("pfam") or [])
    target_kos = set(target.get("ko") or [])
    logic = (target.get("pfam_logic") or "any").lower()   # any | all
    fallback = bool(target.get("blast_fallback"))
    requires_blast_for_confirmation = bool(
        target.get("requires_blast_for_confirmation"))
    pident_min = float(target.get("blast_identity_min")
                       or pident_default)

    # Collect per-protein evidence summaries.
    candidates = []
    proteins = (set(pfam_hits.keys()) | set(custom_hits.keys())
                | set(ko_hits.keys()) | set(blast_hits.keys()))
    for prot in proteins:
        p_hits = pfam_hits.get(prot, {})
        c_hits = custom_hits.get(prot, {})
        k_hits = ko_hits.get(prot, {})
        b_hits_all = blast_hits.get(prot, {}).get(tid, [])
        b_hits = [(a, pi, ev) for (a, pi, ev) in b_hits_all
                  if pi >= pident_min]
        # Sequence-signature check. Precedence: custom HMM > KOfam KO > Pfam.
        # The custom-HMM tier is preferred when a clade model has been built
        # (hardening phase); until then KO is the primary signature, so a
        # `custom_hmm: true` target with no built model still resolves via KO.
        sig_source = None
        sig_evalue = None
        present_pfams: set[str] = set()
        present_kos: set[str] = set()
        if is_custom and tid in c_hits:
            sig_source = "custom-hmm"
            sig_evalue = c_hits[tid]
        if sig_source is None and target_kos:
            present_kos = target_kos & set(k_hits.keys())
            if present_kos:
                sig_source = "ko"
                sig_evalue = min(k_hits[k] for k in present_kos)
        if sig_source is None and req_pfams:
            present_pfams = req_pfams & set(p_hits.keys())
            ok = ((logic == "all" and present_pfams == req_pfams)
                  or (logic == "any" and bool(present_pfams)))
            if ok:
                sig_source = "pfam"
                sig_evalue = min(p_hits[p] for p in present_pfams)
        sig_ok = sig_source is not None
        blast_ok = bool(b_hits)
        if not (sig_ok or blast_ok):
            continue
        # Status.
        # Custom-HMM-confirmed calls bypass the BLAST gate: the clade-calibrated
        # HMM is a stronger discriminator than the BLAST identity floor (e.g.,
        # nxrA HMM at 1979 ≫ calibrated TC 700 ≫ denitrifier band ~530 — see
        # REPORT P3). KO-only and Pfam-only paths still require BLAST confirmation
        # where the target asks for it.
        if (requires_blast_for_confirmation and sig_ok and not blast_ok
                and sig_source != "custom-hmm"):
            status = "disqualified"
        elif sig_ok and (fallback and blast_ok or not fallback):
            status = "confirmed"
        elif sig_ok and fallback and not blast_ok:
            status = "domain-only"
        elif not sig_ok and blast_ok:
            # For custom-HMM targets, the HMM is the authoritative diagnostic.
            # A BLAST-only hit (HMM rejected via TC) is cross-reactivity to a
            # generic homologue, not a real call. For Pfam or Tier-3 targets,
            # narrow-no-IPR remains a meaningful "BLAST says yes, family
            # signature says no" status (e.g. lanmodulin Tier-3 BLAST-only).
            if is_custom and not target_kos:
                continue
            status = "narrow-no-IPR"
        else:
            continue
        evidence_source = sig_source if sig_ok else "blast"
        best_b = min(b_hits, key=lambda x: x[2]) if b_hits else None
        candidates.append({
            "protein_id": prot,
            "status": status,
            "evidence_source": evidence_source,
            "pfam_hits": (sorted(present_kos) if sig_source == "ko"
                          else sorted(present_pfams) if sig_source == "pfam"
                          else [tid] if sig_source == "custom-hmm"
                          else []),
            "pfam_evalue": sig_evalue,
            "blast_acc": best_b[0] if best_b else "",
            "blast_pident": best_b[1] if best_b else "",
            "blast_evalue": best_b[2] if best_b else "",
        })
    if not candidates:
        return None
    # Prefer confirmed > domain-only > narrow-no-IPR > disqualified.
    rank = {"confirmed": 0, "domain-only": 1,
            "narrow-no-IPR": 2, "disqualified": 3}
    best = sorted(candidates,
                  key=lambda c: (rank[c["status"]],
                                 c["pfam_evalue"] or c["blast_evalue"] or 1e9))[0]
    return best["status"], best


# ───────────────────────────── nxr/nar synteny ───────────────────────────────

def parse_prodigal_coords(faa_path: Path | None) -> dict[str, tuple[str, int, int]]:
    """Parse gene coordinates from a Prodigal .faa. Prodigal headers are
    `>{seqid}_{n} # {start} # {end} # {strand} # ID=...`; the downstream protein
    id is `{seqid}_{n}` and the contig is the seqid (id minus the trailing _{n}).
    Returns {protein_id: (contig, start, end)}, or {} for a pre-called proteome
    (no `#`-delimited coordinate header) — synteny resolution is then skipped."""
    out: dict[str, tuple[str, int, int]] = {}
    if not faa_path or not Path(faa_path).exists():
        return out
    with open(faa_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = [p.strip() for p in line[1:].split("#")]
            if len(parts) < 4:
                continue                      # not a Prodigal coordinate header
            protid = parts[0].split()[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            out[protid] = (protid.rsplit("_", 1)[0], start, end)
    return out


def resolve_nxr_synteny(coords: dict[str, tuple[str, int, int]],
                        ko_by_protein: dict[str, dict[str, float]],
                        window: int = 10000) -> dict[str, str]:
    """Disambiguate the K00370/K00371 homology trap (nxrA/narG, nxrB/narH) by
    operon context. A nitrate-reductase alpha (K00370) / beta (K00371) ORF that
    is syntenic with narI (K00374 — the respiratory NAR membrane cytochrome, which
    NXR lacks) is NarG/NarH; one with no adjacent narI is the nitrite-oxidoreductase
    NxrA/NxrB. Nitrobacter NxrA is ~60% identical to NarG and scores in the same
    band, so neither HMM nor BLAST separates them — but the genome context does
    (the narG-narH-narJ-narI operon vs the standalone nxr). Returns {target_id:
    protein_id} for whichever of nxrA/narG/nxrB/narH are present.
    Requires gene coordinates → nucleotide/MAG input only (empty for proteomes)."""
    def narI_within(prot: str, markers: list[str]) -> bool:
        c, s, e = coords[prot]
        lo1, hi1 = min(s, e), max(s, e)
        for m in markers:
            cm, sm, em = coords[m]
            if cm != c:
                continue
            lo2, hi2 = min(sm, em), max(sm, em)
            if max(lo1, lo2) - min(hi1, hi2) <= window:   # gap ≤ window (overlap → <0)
                return True
        return False

    narI = [p for p in coords if "K00374" in ko_by_protein.get(p, {})]
    out: dict[str, str] = {}
    for ko, (operon_tid, nxr_tid) in (("K00370", ("narG", "nxrA")),
                                      ("K00371", ("narH", "nxrB"))):
        for p in (q for q in coords if ko in ko_by_protein.get(q, {})):
            out.setdefault(operon_tid if narI_within(p, narI) else nxr_tid, p)
    return out


# ───────────────────────────── main ──────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", choices=["protein"], default="protein",
                    help="Only protein mode is supported.")
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--hmm", type=Path)
    ap.add_argument("--blast-unstable", type=Path)
    ap.add_argument("--blast-gated", type=Path)
    ap.add_argument("--tc-cutoffs", type=Path,
                    default=Path("resources/hmm/tc_cutoffs.tsv"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pident-default", type=float, default=30.0)
    ap.add_argument("--qcov-default", type=float, default=50.0)
    ap.add_argument("--gene-coords", type=Path,
                    help="proteome .faa; if it carries Prodigal coordinate headers "
                         "(nucleotide/MAG input), resolve the nxrA/narG + nxrB/narH "
                         "trap by operon synteny. Ignored for pre-called proteomes.")
    args = ap.parse_args()

    with open(args.targets) as fh:
        targets = yaml.safe_load(fh)["targets"]

    if args.mode == "protein":
        tc = load_tc_cutoffs(args.tc_cutoffs)
        # Per-target hmm_min_qcov override from targets.yaml — defends small
        # domain HMMs (e.g. 104-aa hcnA) from hitting domains of large
        # multidomain proteins.
        qcov_cutoffs = {t["id"]: float(t["hmm_min_qcov"])
                        for t in targets if t.get("hmm_min_qcov") is not None}
        hmm_parsed = parse_hmmscan_domtbl(args.hmm, tc, qcov_cutoffs) if args.hmm else {}
        blast_hits = {}
        for p in (args.blast_unstable, args.blast_gated):
            if p and p.exists():
                for prot, td in parse_blast_tsv(p).items():
                    blast_hits.setdefault(prot, {}).update(td)
        # Split per signature tier: pfam_hits[prot]={pfam_id: evalue};
        # custom_hits[prot]={tid: evalue}; ko_hits[prot]={KO: evalue}.
        pfam_hits_flat = {p: d["pfam"] for p, d in hmm_parsed.items()}
        custom_hits_flat = {p: d["custom"] for p, d in hmm_parsed.items()}
        ko_hits_flat = {p: d["ko"] for p, d in hmm_parsed.items()}

        rows = []
        for t in targets:
            res = evaluate_target(t, pfam_hits_flat, custom_hits_flat,
                                  ko_hits_flat, blast_hits,
                                  args.pident_default, args.qcov_default)
            if res is None:
                rows.append({
                    "target_id": t["id"], "name": t["name"],
                    "category": t["category"], "complex": t.get("complex", ""),
                    "status": "absent", "evidence_source": "",
                    "protein_id": "", "pfam_hits": "",
                    "best_pfam_evalue": "", "blast_acc": "",
                    "blast_pident": "", "blast_evalue": "",
                })
                continue
            status, ev = res
            rows.append({
                "target_id": t["id"], "name": t["name"],
                "category": t["category"], "complex": t.get("complex", ""),
                "status": status,
                "evidence_source": ev.get("evidence_source", ""),
                "protein_id": ev["protein_id"],
                "pfam_hits": ",".join(ev["pfam_hits"]),
                "best_pfam_evalue": ev["pfam_evalue"],
                "blast_acc": ev["blast_acc"],
                "blast_pident": ev["blast_pident"],
                "blast_evalue": ev["blast_evalue"],
            })

        # Operon-synteny resolution of the nxrA/narG + nxrB/narH trap. Only fires
        # when --gene-coords carries Prodigal coordinates (nucleotide/MAG input);
        # for pre-called proteomes coords are {} and the calls above stand.
        coords = parse_prodigal_coords(args.gene_coords)
        if coords:
            syn = resolve_nxr_synteny(coords, ko_hits_flat)
            by_id = {r["target_id"]: r for r in rows}
            for tid in ("nxrA", "narG", "nxrB", "narH"):
                r = by_id.get(tid)
                if r is None:
                    continue
                if tid in syn:
                    r.update(status="confirmed", evidence_source="synteny",
                             protein_id=syn[tid], pfam_hits=tid,
                             best_pfam_evalue="", blast_acc="",
                             blast_pident="", blast_evalue="")
                elif not (r.get("evidence_source") == "custom-hmm"
                          and r["status"] != "absent"):
                    # Synteny gave no verdict for this tid. Preserve a present
                    # custom-HMM call (confirmed OR domain-only) — e.g. comammox
                    # NxrA (≫ TC 700, BLAST-confirmed) and Nitrospina NxrA (≫ TC but
                    # divergent periplasmic NXR fails the BLAST gate → domain-only):
                    # the clade HMM is a stronger discriminator than synteny and must
                    # not be wiped when synteny is silent. Otherwise clear the
                    # trap-ambiguous KO/Pfam-only call (the conservative non-call for
                    # the sequence-inseparable Nitrobacter NxrA↔NarG case).
                    r.update(status="absent", evidence_source="", protein_id="",
                             pfam_hits="", best_pfam_evalue="", blast_acc="",
                             blast_pident="", blast_evalue="")
            print(f"[apply_rules] {args.sample}: synteny-resolved nxr/nar trap: "
                  f"{syn or '{}'}", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    n_status = defaultdict(int)
    for r in rows:
        n_status[r["status"]] += 1
    print(f"[apply_rules] {args.sample}: " +
          ", ".join(f"{k}={v}" for k, v in sorted(n_status.items())),
          file=sys.stderr)


if __name__ == "__main__":
    main()
