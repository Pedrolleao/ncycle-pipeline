#!/usr/bin/env python3
"""
adapters.py — load each tool's calls into ONE shared vocabulary, plus ground truth.

Every loader returns subunit-level predictions: {(genome, target_id): bool_present},
restricted to the training panel (hold-outs excluded, as in compare_kofam.py). The
benchmark engine collapses to the `step` resolution itself (STEP_DIAG below), so the
collapse rule is applied identically to every tool and to the ground truth.

Working now: load_truth, load_ncycle (our pipeline), load_kofam (raw KO baseline).
Scaffolded: load_metabolic / load_ncycdb / load_dram — each parses a normalized TSV
exported from that tool's native output and maps its identifiers to our targets. They
return {} when given no input path, so the engine simply skips a tool until its
output is provided. Fill the mapping tables + the exact tool version/command when run.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GT = ROOT / "validation" / "ground_truth.tsv"
MATRIX = ROOT / "ncycle_results" / "ncycle_matrix.tsv"
TARGETS = ROOT / "config" / "targets.yaml"
KO_LIST = ROOT / "resources" / ".cache" / "ko_list"
RESULTS = ROOT / "ncycle_results"

HOLDOUT_TAG = "holdout_v3"

# Homology-trap targets — the confirmatory-claim subset (subunit resolution).
TRAP = {"nxrA", "nxrB", "narG", "narH", "amoA", "amoB", "amoC", "napA",
        "nirK", "nirB", "nirD", "nirA", "amoA_archaeal"}

# ── step resolution ──────────────────────────────────────────────────────────
# STEP_MEMBERS: every subunit that belongs to a step (used to map a coarse tool's
# step-level call DOWN to subunits). STEP_DIAG: the catalytic/diagnostic marker(s)
# whose presence defines the step (used to collapse subunit calls + truth UP to steps).
STEP_MEMBERS: dict[str, list[str]] = {
    "N2_fixation":              ["nifH","nifD","nifK","nifE","nifN","nifB","vnfH","vnfD","anfG"],
    "ammonia_oxidation":        ["amoA","amoA_archaeal","amoB","amoC","hao"],
    "nitrite_oxidation":        ["nxrA","nxrB"],
    "dissim_nitrate_reduction": ["narG","narH","narI","napA","napB"],
    "NO2_to_NO":                ["nirK","nirS"],
    "NO_reduction":             ["norB","norC","norZ"],
    "N2O_reduction":            ["nosZ"],
    "DNRA":                     ["nrfA","nrfH","nirB","nirD"],
    "anammox":                  ["hzsA","hzsB","hzsC","hdh"],
    "assim_nitrate_reduction":  ["narB","nasA","nasB"],
    "assim_nitrite_reduction":  ["nirA","nasD","NR"],
    "ammonia_assimilation":     ["glnA","gltB","gltD","gdhA"],
    "ureolysis":                ["ureA","ureB","ureC","ureG"],
}
STEP_DIAG: dict[str, list[str]] = {
    "N2_fixation":              ["nifH","nifD"],
    "ammonia_oxidation":        ["amoA","amoA_archaeal"],
    "nitrite_oxidation":        ["nxrA"],
    "dissim_nitrate_reduction": ["narG","napA"],
    "NO2_to_NO":                ["nirK","nirS"],
    "NO_reduction":             ["norB","norZ"],
    "N2O_reduction":            ["nosZ"],
    "DNRA":                     ["nrfA","nirB"],
    "anammox":                  ["hzsA","hdh"],
    "assim_nitrate_reduction":  ["narB","nasA"],
    "assim_nitrite_reduction":  ["nirA","nasD","NR"],
    "ammonia_assimilation":     ["glnA","gdhA"],
    "ureolysis":                ["ureC"],
}
TARGET_STEP = {t: step for step, members in STEP_MEMBERS.items() for t in members}


# ── shared helpers ───────────────────────────────────────────────────────────

def _cfg_targets() -> list[dict]:
    import yaml
    return yaml.safe_load(open(TARGETS))["targets"]


def target_kos() -> dict[str, list[str]]:
    return {t["id"]: (t.get("ko") or []) for t in _cfg_targets()}


def all_target_ids() -> list[str]:
    return [t["id"] for t in _cfg_targets()]


def load_truth() -> tuple[dict[tuple[str, str], str], set[str]]:
    """Return ({(genome,target): 'present'|'absent'}, training_genomes), hold-outs dropped."""
    truth, genomes = {}, set()
    with open(GT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("source") == HOLDOUT_TAG:
                continue
            truth[(r["genome"], r["target"])] = r["expected"]
            genomes.add(r["genome"])
    return truth, genomes


# ── tool loaders (subunit-level predictions) ─────────────────────────────────

def load_ncycle(genomes: set[str]) -> dict[tuple[str, str], bool]:
    """Our pipeline: matrix status codes 1/2 = present."""
    pred: dict[tuple[str, str], bool] = {}
    with open(MATRIX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["sample"]
            if g not in genomes:
                continue
            for col, val in r.items():
                if col.startswith("target__"):
                    try:
                        pred[(g, col[len("target__"):])] = int(val) in (1, 2)
                    except (TypeError, ValueError):
                        pred[(g, col[len("target__"):])] = False
    return pred


def load_kofam(genomes: set[str]) -> dict[tuple[str, str], bool]:
    """Raw KofamScan baseline: target present iff ANY of its KOs clears the stock
    ko_list threshold (no gating, no custom HMMs; shared KOs → all their targets)."""
    ko_thr: dict[str, float] = {}
    with open(KO_LIST) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[1] not in ("", "-"):
                try:
                    ko_thr[p[0]] = float(p[1])
                except ValueError:
                    pass
    tkos = target_kos()
    pred: dict[tuple[str, str], bool] = {}
    for g in genomes:
        tbl = RESULTS / g / "hmm" / f"{g}.hmmscan.tsv"
        present_kos = set()
        if tbl.exists():
            with open(tbl) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    c = line.split()
                    if len(c) < 8:
                        continue
                    # KO name + full-sequence score column depend on the HMMER layout:
                    # hmmsearch --domtblout (current pipeline) → query=col3, score=col7;
                    # legacy hmmscan --tblout → target=col0, score=col5. Detect by which
                    # column holds a known KO so the raw-KO baseline survives either format.
                    if c[3] in ko_thr:
                        ko, score = c[3], c[7]
                    elif c[0] in ko_thr:
                        ko, score = c[0], c[5]
                    else:
                        continue
                    try:
                        if float(score) >= ko_thr[ko]:
                            present_kos.add(ko)
                    except ValueError:
                        pass
        for tid, kos in tkos.items():
            pred[(g, tid)] = any(k in present_kos for k in kos)
    return pred


# ── external comparator adapters ──
# RUN PROVENANCE — all on the 31 training proteomes (comparators/{metabolic,dram}_in/, copied from
# ../test_panel), default settings; tool installs reused from Sulfur_Cycle/comparators.
#   METABOLIC v4.0 (env METABOLIC_v4.0) — RUN 2026-06-09:
#     `perl METABOLIC-G.pl -t 8 -in comparators/metabolic_in -kofam-db full -o comparators/metabolic_out`.
#     metabolic.tsv = (genome, ko): a KO is present iff METABOLIC's per-genome KEGG result
#     (metabolic_out/KEGG_identifier_result/<genome>.result.txt) lists a protein hit for it. KO→target
#     via the same map as raw kofam, so METABOLIC's trap behavior is governed by the same shared-KO
#     ambiguity (it under-calls due to stricter KOfam m-cutoff). Built by comparators/build_metabolic_tsv.py.
#     Its bundled KOfam lacks K10534 (NR) + K17877 (nasD) — 2 KO-version coverage gaps.
#   DRAM v1.4.6 (env DRAM14; KOfam-only config — Pfam/dbCAN/MEROPS = None, KO-irrelevant here) —
#     RUN 2026-06-09: `DRAM.py annotate_genes -i 'comparators/dram_in/*.faa' -o comparators/dram_out
#     --threads 8` then `DRAM.py distill`. dram.tsv = (genome, module, present): a KEGG module
#     (DRAM_STEP_MAP below) is present iff DRAM detected ≥1 of its member KOs (from the
#     genome_summary_form) in that genome's annotations.tsv. Module→subunit expansion makes DRAM
#     coarsely over-call the traps; it has no module for glnA/gltBD/gdhA or catalytic ureC (gaps).
#     Built by comparators/build_dram_tsv.py. (dram_in headers genome-prefixed for unique gene ids.)
#   NCycDB (Tu et al. 2019, github.com/qichao1984/NCyc) — RUN 2026-06-09. NCyc_100.faa (219,146
#     seqs @100% id) + id2gene.map (68 gene families). DIAMOND 2.1.9, faithful to NCycProfiler.PL
#     protein default: `diamond makedb --in NCyc_100.faa --db NCyc_100` then per training proteome
#     `diamond blastp -k 1 -e 0.0001 -d NCyc_100 -q <faa>`. Best-hit subject -> family via
#     id2gene.map; ncycdb.tsv = (sample, family, count). Built by ../../comparators/build_ncycdb_tsv.py
#     (DB + per-genome hits under comparators/NCyc/ + comparators/ncycdb_out/).
#     amoA_A/amoA_B were verified empirically clade-specific on the panel: amoA_A is hit ONLY by the
#     two AOA archaea (Nmaritimus, Nviennensis) and amoA_B only by bacterial AOB/comammox/gamma-AOB
#     (+ a spurious Paeruginosa best-hit) -> amoA_A=amoA_archaeal, amoA_B=amoA. NCyc has no gamma-AOB
#     amoA family (Noceani's gamma-amoA falls into amoA_B), no nrfH, no qNor norZ, no nosZ clade-II,
#     and no nifE/N/B or vnf accessory families -> those ncycle targets are NCyc coverage gaps (FN).

def _read_norm_tsv(path: Path | None) -> list[dict]:
    """Read a normalized 'present-call' TSV exported from a tool's native output."""
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_metabolic(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    """METABOLIC (Zhou et al. 2022). METABOLIC reports KO/function presence per genome;
    export it to a normalized TSV with columns: genome, ko (one present KO per row).
    A KO present → all its ncycle targets present (same KO→target map as KofamScan).
    Returns {} if no path. Verify the export against METABOLIC_result.xlsx when run."""
    rows = _read_norm_tsv(path)
    if not rows:
        return {}
    tkos = target_kos()
    present_kos: dict[str, set[str]] = {}
    for r in rows:
        present_kos.setdefault(r["genome"], set()).add(r["ko"])
    return {(g, tid): any(k in present_kos.get(g, set()) for k in kos)
            for g in genomes for tid, kos in tkos.items()}


# NCyc gene-family name → ncycle target id. COMPLETE over NCyc's 68 families (those with no
# corresponding ncycle target are intentionally omitted → NCyc predicts them absent; those are
# either non-targets, e.g. pmoA/ansB/hcp, or genuine NCyc coverage gaps documented in the
# provenance block above). Mapping fixed from NCyc's family semantics + the empirical amoA_A/B
# clade check; applied identically to all genomes (prereg: vocabulary is a fixed, reviewable table).
NCYC_MAP: dict[str, str] = {
    # ammonia oxidation — NCyc resolves clade: _A = archaeal (AOA), _B = bacterial (AOB/comammox)
    "amoA_A": "amoA_archaeal", "amoA_B": "amoA",
    "amoB_A": "amoB", "amoB_B": "amoB", "amoC_A": "amoC", "amoC_B": "amoC", "hao": "hao",
    # nitrite oxidation (NCyc keeps nxr separate from nar, like ncycle)
    "nxrA": "nxrA", "nxrB": "nxrB",
    # dissimilatory nitrate reduction
    "narG": "narG", "narH": "narH", "narI": "narI", "napA": "napA", "napB": "napB",
    # NO2 -> NO ; NO reduction ; N2O reduction
    "nirK": "nirK", "nirS": "nirS", "norB": "norB", "norC": "norC", "nosZ": "nosZ",
    # DNRA (NCyc has nrfA + nirBD; nrfH is a coverage gap)
    "nrfA": "nrfA", "nirB": "nirB", "nirD": "nirD",
    # N2 fixation (NCyc has nifHDK + anfG only; nifE/N/B + vnf are coverage gaps)
    "nifH": "nifH", "nifD": "nifD", "nifK": "nifK", "anfG": "anfG",
    # assimilatory nitrate/nitrite reduction
    "narB": "narB", "nasA": "nasA", "nasB": "nasB", "nirA": "nirA", "NR": "NR",
    # ammonia assimilation (gs/gdh families keyed by KO; map only our gltB/gltD/gdhA KOs)
    "glnA": "glnA", "gs_K00265": "gltB", "gs_K00266": "gltD",
    "gdh_K00261": "gdhA", "gdh_K00262": "gdhA",
    # anammox
    "hzsA": "hzsA", "hzsB": "hzsB", "hzsC": "hzsC", "hdh": "hdh",
    # ureolysis
    "ureA": "ureA", "ureB": "ureB", "ureC": "ureC",
}


def load_ncycdb(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    """NCycDB (Tu et al. 2019) + DIAMOND. Export a normalized TSV: columns
    sample, family, count. family→target via NCYC_MAP; present iff count>0.
    Returns {} if no path."""
    rows = _read_norm_tsv(path)
    if not rows:
        return {}
    hit: dict[str, set[str]] = {}
    for r in rows:
        tid = NCYC_MAP.get(r["family"])
        if tid and float(r.get("count", 1) or 0) > 0:
            hit.setdefault(r["genome"] if "genome" in r else r["sample"], set()).add(tid)
    return {(g, tid): tid in hit.get(g, set())
            for g in genomes for tid in all_target_ids()}


# DRAM distillate MODULE label → ncycle STEP(s). DRAM reports at module granularity; the
# engine expands a present module to ALL member subunits of EVERY step it covers (symmetric
# credit/blame) — this is the honest representation of what DRAM tells a user, and the source
# of its coarse over-calling on the homology traps.
#
# Labels are DRAM 1.4.6's ACTUAL `genome_summary_form` 'Energy'/Nitrogen module names (verified
# 2026-06-07 against DRAM_data/forms/genome_summary_form.20260527.tsv). DRAM lumps the entire
# denitrification chain into one module (nar→nir→nor→nos share it) and cannot separate
# nitrite-oxidation (nxr) from nitrate-reduction (nar) — they share K00370/K00371 — so the
# comammox/denitrification/DNRA modules all carry those KOs. DRAM has NO module for the catalytic
# urease (ureC = K01428-30; its "Urea cycle" module is the arginine/urea-cycle, K00611/K01476/…)
# nor for ammonia assimilation (glnA/gltBD/gdhA) — so ureolysis + ammonia_assimilation are scored
# absent for DRAM (a genuine DRAM coverage limitation, not a mapping gap).
DRAM_STEP_MAP: dict[str, list[str]] = {
    "nitrogen fixation, nitrogen => ammonia": ["N2_fixation"],
    "nitrification, ammonia => nitrite": ["ammonia_oxidation"],
    "complete nitrification, comammox, ammonia => nitrite => nitrate":
        ["ammonia_oxidation", "nitrite_oxidation"],
    "denitrification, nitrate => nitrogen":
        ["dissim_nitrate_reduction", "NO2_to_NO", "NO_reduction", "N2O_reduction"],
    "dissimilatory nitrate reduction, nitrate => ammonia":
        ["dissim_nitrate_reduction", "DNRA"],
    "assimilatory nitrate reduction, nitrate => ammonia":
        ["assim_nitrate_reduction", "assim_nitrite_reduction"],
    "nitrate assimilation": ["assim_nitrate_reduction"],
    "nitrite + ammonia => nitrogen": ["anammox"],
}


def load_dram(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    """DRAM (Shaffer et al. 2020). Reads a normalized TSV: columns genome, function, present
    (function = DRAM module label; present iff DRAM detected ≥1 of the module's genes in that
    genome). A present module is expanded to all STEP_MEMBERS subunits of EVERY step it maps to.
    Returns {} if no path."""
    rows = _read_norm_tsv(path)
    if not rows:
        return {}
    step_present: dict[str, set[str]] = {}
    for r in rows:
        steps = DRAM_STEP_MAP.get(r["function"].strip().lower(), [])
        if steps and str(r.get("present", "")).lower() in ("1", "true", "yes", "present"):
            step_present.setdefault(r["genome"], set()).update(steps)
    pred: dict[tuple[str, str], bool] = {}
    for g in genomes:
        present_subunits = {t for s in step_present.get(g, set()) for t in STEP_MEMBERS[s]}
        for tid in all_target_ids():
            pred[(g, tid)] = tid in present_subunits
    return pred


# Registry: name → (loader, needs_path). The engine calls each; empty result = skipped.
LOADERS = {
    "ncycle":   (load_ncycle,   False),
    "kofam":    (load_kofam,    False),
    "metabolic": (load_metabolic, True),
    "ncycdb":   (load_ncycdb,   True),
    "dram":     (load_dram,     True),
}
