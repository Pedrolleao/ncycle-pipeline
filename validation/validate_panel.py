#!/usr/bin/env python3
"""
validate_panel.py — per-genome QC gate for the reference panel.

Added 2026-05-30 (ROADMAP P5.0.6) after the 2026-05-30 audit found
`test_panel/Nwinogradskyi_Nb255.faa` was a *Burkholderia thailandensis* proteome
(5258/5607 headers `[Burkholderia thailandensis]`-tagged), not Nitrobacter
winogradskyi. The contaminated file invalidated 3 Nb-255 GT corrections
(v6 nosZ, v9 norZ qNor, and stress-tested v4 narG); see
validation/REPORT.md § Audit 2026-05-30 → Post-audit response (v11).

Gate logic
----------
For each .faa in the panel:

  1. Sample the first N protein headers (default 50; or all if smaller).
  2. For each header, extract an organism string (UniProt `OS=…` or NCBI
     `[Taxon]`) and classify against the expected taxa:
        - `match`     — organism string contains an expected taxon token
        - `wrong`     — organism string names a specific binomial *other than*
                        the expected one (the Burkholderia signal)
        - `agnostic`  — header is rank-elevated (`[Bacteria]`, `[Pseudomonadaceae]`,
                        `[Gammaproteobacteria]`, MULTISPECIES …) or carries no
                        taxon info. Cannot confirm or refute the expectation;
                        does not count toward the wrong-genus penalty.
  3. **Fail** the genome if `wrong_count / max(match_count + wrong_count, 1)
     > WRONG_GENUS_TOLERANCE` (default 0.05 — i.e. >5% of identifiable-
     species headers name an unexpected genus).
  4. For prodigal-derived proteomes (no organism tags in headers, only contig
     IDs like `NC_007406.1_N`), use the contig-prefix path instead: require
     ≥95% of sampled headers to start with one of the expected contig
     prefixes.

The thresholds are deliberately permissive on `agnostic` headers so we don't
penalize MULTISPECIES proteomes with heavy rank-elevation, while still
catching a Burkholderia-style swap where the *majority of organism-tagged
headers explicitly name a different genus*.

Standalone usage
----------------
    python validation/validate_panel.py            # default panel dir
    python validation/validate_panel.py --panel ../test_panel
    python validation/validate_panel.py --strict   # fail on any agnostic-only genome

Wired into the Makefile via the `validate-panel` target; `regression` depends on it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PANEL = (ROOT / ".." / "test_panel").resolve()

SAMPLE_SIZE = 50
WRONG_GENUS_TOLERANCE = 0.05   # >5% wrong-genus → fail
CONTIG_MATCH_FLOOR = 0.95      # ≥95% contig-prefix → pass on prodigal-mode files

# Per-genome expectation table.
#
# `taxa`: lowercase tokens. A header organism string is `match` if it contains
#   ANY of these tokens as a whole word. Include genus + species; add reclassified
#   names (e.g. Desulfovibrio→Nitratidesulfovibrio) as accepted synonyms.
# `contig_prefix`: optional list. If set AND the file has no organism-tagged
#   headers (prodigal-derived), require headers to start with one of these.
# `note`: free-text caveat surfaced in the report (e.g. file vs. content species
#   mismatch that has been verified safe).
EXPECTED: dict[str, dict] = {
    "Ahydrophila_ATCC7966":      {"taxa": ["aeromonas", "hydrophila"]},
    "Aterreus_NIH2624":          {"taxa": ["aspergillus", "terreus"]},
    "Avinelandii_DJ":            {"taxa": ["azotobacter", "vinelandii"]},
    "Bdiazoefficiens_USDA110":   {"taxa": ["bradyrhizobium", "diazoefficiens",
                                            "japonicum"]},  # historical synonym
    "Bsinica_JPN1":              {"taxa": ["brocadia", "sinica"]},
    "Bsubtilis_168":             {"taxa": ["bacillus", "subtilis"]},
    "Cnecator_H16":              {"taxa": ["cupriavidus", "necator",
                                            "ralstonia", "wautersia"]},  # historical synonyms
    "Dvulgaris_Hildenborough":   {"taxa": ["desulfovibrio", "nitratidesulfovibrio",
                                            "vulgaris"]},  # NCBI reclassified
    "Ecoli_K12_MG1655":          {"taxa": ["escherichia", "coli"]},
    "Hpylori_26695":             {"taxa": ["helicobacter", "pylori"]},
    "Kstuttgartiensis":          {"taxa": ["kuenenia", "stuttgartiensis"]},
    "Lacidophilus_4356":         {"taxa": ["lactobacillus", "acidophilus"]},
    "Mfumariolicum_SolV":        {"taxa": ["methylacidiphilum", "fumariolicum"]},
    "Neuropaea_ATCC19718":       {"taxa": ["nitrosomonas", "europaea"]},
    "Ngracilis_3211":            {"taxa": ["nitrospina", "gracilis"]},
    "Ninopinata_comammox":       {"taxa": ["nitrospira", "inopinata"]},
    "Njaponica_NJ11":            {"taxa": ["nitrospira", "japonica"]},
    "Nmaritimus_SCM1":           {"taxa": ["nitrosopumilus", "maritimus"]},
    "Nmultiformis_ATCC25196":    {"taxa": ["nitrosospira", "multiformis"]},
    "Nnitrosa_comammox":         {"taxa": ["nitrospira", "nitrosa", "nitrificans"],
                                   "note": "panel file is filename-tagged 'nitrosa' "
                                           "but content headers cite 'Candidatus "
                                           "Nitrospira nitrificans' — same comammox "
                                           "genus, different species. Accepted."},
    "Noceani_ATCC19707":         {"taxa": ["nitrosococcus", "oceani"]},
    "Npcc7120":                  {"taxa": ["nostoc", "anabaena"]},  # genus has been split
    "Nviennensis_EN76":          {"taxa": ["nitrososphaera", "viennensis"]},
    "Nwinogradskyi_Nb255":       {"taxa": ["nitrobacter", "winogradskyi"],
                                   "contig_prefix": ["NC_007406.1"],
                                   "note": "prodigal-derived from test_synteny/"
                                           "Nwinogradskyi.fna (NC_007406.1); headers "
                                           "carry no organism tags, contig-prefix "
                                           "fallback used."},
    "Paeruginosa_PAO1":          {"taxa": ["pseudomonas", "aeruginosa"]},
    "Pdenitrificans_PD1222":     {"taxa": ["paracoccus", "denitrificans"]},
    "Rpalustris_CGA009":         {"taxa": ["rhodopseudomonas", "palustris"]},
    "Scerevisiae_S288C":         {"taxa": ["saccharomyces", "cerevisiae"]},
    "Smeliloti_1021":            {"taxa": ["sinorhizobium", "ensifer", "meliloti"]},
    "Soneidensis_MR1":           {"taxa": ["shewanella", "oneidensis"]},
    "Spneumoniae_ref":           {"taxa": ["streptococcus", "pneumoniae"]},
    "Synechocystis_PCC6803":     {"taxa": ["synechocystis"]},
    "Wsuccinogenes_DSM1740":     {"taxa": ["wolinella", "succinogenes"]},
    # P5.3.5 (2026-05-30): 6 new hold-out genomes spanning 5 phyla.
    "Adehalogenans_2CP1":        {"taxa": ["anaeromyxobacter", "dehalogenans"]},
    "Nhalophilus_Nc4":           {"taxa": ["nitrosococcus", "halophilus"]},
    "Nhollandica_Lb":            {"taxa": ["nitrolancea", "hollandica"]},
    "Sbrodae":                   {"taxa": ["scalindua", "brodae"]},
    "Sstutzeri_F2a":             {"taxa": ["stutzerimonas", "pseudomonas", "stutzeri"],
                                   "note": "NCBI reclassified Pseudomonas stutzeri → "
                                           "Stutzerimonas stutzeri; accept either genus."},
    "Mcapsulatus_Bath":          {"taxa": ["methylococcus", "capsulatus"]},
}

# Higher-rank suffixes — markers of agnostic (non-binomial) organism strings.
# A header organism whose first word ends in one of these is rank-elevated and
# cannot be used to confirm OR refute the genus expectation.
HIGHER_RANK_SUFFIXES = (
    "aceae",        # families: Pseudomonadaceae, Nitrosomonadaceae
    "ales",         # orders: Burkholderiales
    "ineae",        # suborders
    "idae",         # subclasses
    "phyceae",      # algal classes: Cyanophyceae
    "phyta",        # plant phyla
    "mycota",       # fungal phyla
    "mycetes",      # fungal classes
    "bacteriaceae",
)

# Specific high-rank token names (single-word taxa above genus that are not
# matched by suffix rules — domain, kingdom, phylum, class).
HIGHER_RANK_NAMES = {
    "bacteria", "archaea", "eukaryota", "viruses",
    "proteobacteria", "actinobacteriota", "actinobacteria", "firmicutes",
    "bacillota", "bacteroidota", "bacteroidetes",
    "cyanobacteria", "cyanobacteriota",
    "pseudomonadota", "pseudomonadati", "spirochaetes", "spirochaetota",
    "verrucomicrobiota", "verrucomicrobia",
    "planctomycetota", "planctomycetes",
    "thermotogota", "chloroflexota", "chloroflexi", "fusobacteriota",
    "alphaproteobacteria", "betaproteobacteria", "gammaproteobacteria",
    "deltaproteobacteria", "epsilonproteobacteria", "zetaproteobacteria",
    "campylobacterota",
    "fungi", "metazoa", "viridiplantae",
}

# Tokens that should not be treated as genus names when they happen to be
# the first word of an organism string ("Candidatus" is a status not a name).
NON_TAXON_PREFIX = {"candidatus"}

UNIPROT_OS_RE = re.compile(r"OS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)")
NCBI_BRACKET_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


def extract_organism(header: str) -> str | None:
    """Pull an organism string from a FASTA defline.

    Tries UniProt-style `OS=...` first, then NCBI-style trailing `[Taxon]`.
    Returns None if neither pattern matches (e.g. prodigal headers).
    """
    m = UNIPROT_OS_RE.search(header)
    if m:
        return m.group(1).strip()
    m = NCBI_BRACKET_RE.search(header)
    if m:
        return m.group(1).strip()
    return None


def classify_organism(organism: str | None, expected_taxa: list[str]) -> str:
    """Return one of: 'match', 'wrong', 'agnostic'."""
    if not organism:
        return "agnostic"
    org_lc = organism.lower()
    # Drop leading 'Candidatus ' — it's a taxonomic status, not a name.
    tokens = [t for t in org_lc.split() if t not in NON_TAXON_PREFIX]
    if not tokens:
        return "agnostic"
    # Whole-word match on any expected taxon.
    for tok in expected_taxa:
        if re.search(r"\b" + re.escape(tok) + r"\b", " ".join(tokens)):
            return "match"
    # Rank-elevation check on the first informative token.
    first = tokens[0]
    if first in HIGHER_RANK_NAMES:
        return "agnostic"
    if any(first.endswith(suf) for suf in HIGHER_RANK_SUFFIXES):
        return "agnostic"
    # Single-token genus or two-token binomial that does NOT match → wrong.
    return "wrong"


def sample_headers(faa: Path, n: int) -> list[str]:
    """Read up to n FASTA defline lines from a .faa file."""
    out: list[str] = []
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line.rstrip())
                if len(out) >= n:
                    break
    return out


def check_genome(faa: Path, expectation: dict, n: int) -> dict:
    """Run the QC checks for a single panel file.

    Returns a dict with: name, status (pass|fail|unknown), match/wrong/agnostic
    counts, contig-prefix metrics if applicable, and a short reason.
    """
    headers = sample_headers(faa, n)
    if not headers:
        return {"name": faa.name, "status": "fail", "reason": "empty file (no headers)",
                "match": 0, "wrong": 0, "agnostic": 0, "n_sampled": 0}

    organisms = [extract_organism(h) for h in headers]
    has_organism = sum(1 for o in organisms if o)

    prefixes = expectation.get("contig_prefix") or []
    # Prodigal-mode path: no organism tags AND we have a contig-prefix expectation.
    if has_organism == 0 and prefixes:
        prefix_hits = sum(
            1 for h in headers
            if any(h.lstrip(">").startswith(p) for p in prefixes)
        )
        prefix_frac = prefix_hits / len(headers)
        if prefix_frac >= CONTIG_MATCH_FLOOR:
            status = "pass"
            reason = (f"contig-prefix match {prefix_hits}/{len(headers)} "
                      f"({prefix_frac:.0%}); expected prefix(es): "
                      f"{', '.join(prefixes)}")
        else:
            status = "fail"
            reason = (f"contig-prefix match {prefix_hits}/{len(headers)} "
                      f"({prefix_frac:.0%}) below {CONTIG_MATCH_FLOOR:.0%} floor")
        return {"name": faa.name, "status": status, "reason": reason,
                "match": prefix_hits, "wrong": len(headers) - prefix_hits,
                "agnostic": 0, "n_sampled": len(headers),
                "mode": "contig-prefix"}

    counts = {"match": 0, "wrong": 0, "agnostic": 0}
    wrong_examples: list[str] = []
    for org in organisms:
        cls = classify_organism(org, expectation["taxa"])
        counts[cls] += 1
        if cls == "wrong" and len(wrong_examples) < 3:
            wrong_examples.append(org or "")

    resolution_bearing = counts["match"] + counts["wrong"]
    if resolution_bearing == 0:
        # All headers were rank-elevated/agnostic — can't confirm OR refute.
        # Treat as 'unknown' (pass with a warning); the panel curator should
        # add a contig_prefix expectation or accept this limitation.
        return {"name": faa.name, "status": "unknown", "reason":
                f"all {len(headers)} sampled headers were rank-elevated or "
                f"untagged; cannot confirm genus from this sample",
                **counts, "n_sampled": len(headers), "mode": "organism-tag"}

    wrong_frac = counts["wrong"] / resolution_bearing
    if wrong_frac > WRONG_GENUS_TOLERANCE:
        status = "fail"
        reason = (f"{counts['wrong']}/{resolution_bearing} resolution-bearing "
                  f"headers name a non-expected genus ({wrong_frac:.0%} > "
                  f"{WRONG_GENUS_TOLERANCE:.0%} tolerance). "
                  f"Examples: {wrong_examples}")
    else:
        status = "pass"
        reason = (f"{counts['match']} match, {counts['wrong']} wrong, "
                  f"{counts['agnostic']} agnostic over {len(headers)} sampled "
                  f"({counts['match']}/{resolution_bearing} resolution-bearing match)")
    return {"name": faa.name, "status": status, "reason": reason,
            **counts, "n_sampled": len(headers), "mode": "organism-tag"}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--panel", type=Path, default=DEFAULT_PANEL,
                   help=f"Panel directory (default {DEFAULT_PANEL})")
    p.add_argument("--sample-size", type=int, default=SAMPLE_SIZE,
                   help=f"Headers to sample per genome (default {SAMPLE_SIZE})")
    p.add_argument("--strict", action="store_true",
                   help="Treat 'unknown' (all-agnostic) genomes as failure")
    args = p.parse_args()

    if not args.panel.is_dir():
        sys.exit(f"error: panel dir not found: {args.panel}")

    files = sorted(args.panel.glob("*.faa"))
    if not files:
        sys.exit(f"error: no .faa files in {args.panel}")

    print(f"[validate_panel] panel: {args.panel}  ({len(files)} .faa files; "
          f"sample size {args.sample_size}; wrong-genus tolerance "
          f"{WRONG_GENUS_TOLERANCE:.0%})")
    print()

    results: list[dict] = []
    untracked: list[str] = []
    for f in files:
        stem = f.stem
        if stem not in EXPECTED:
            untracked.append(f.name)
            continue
        results.append(check_genome(f, EXPECTED[stem], args.sample_size))

    # Print per-genome table.
    n_pass = n_fail = n_unknown = 0
    for r in results:
        marker = {"pass": "  PASS", "fail": "✗ FAIL", "unknown": "? WARN"}[r["status"]]
        print(f"  {marker}  {r['name']:<42s}  {r['reason']}")
        if r["status"] == "pass":   n_pass += 1
        elif r["status"] == "fail": n_fail += 1
        else:                       n_unknown += 1

    if untracked:
        print()
        print(f"[validate_panel] {len(untracked)} file(s) not in EXPECTED map "
              f"(add an entry to validate_panel.py to include in QC):")
        for name in untracked:
            print(f"    {name}")

    print()
    summary = (f"[validate_panel] {n_pass} pass / {n_fail} fail / {n_unknown} unknown "
               f"({len(results)} checked, {len(untracked)} untracked)")
    print(summary)

    if n_fail > 0:
        return 1
    if args.strict and n_unknown > 0:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
