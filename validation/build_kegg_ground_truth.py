#!/usr/bin/env python3
"""
build_kegg_ground_truth.py — KEGG-derived N-cycle ground truth (kegg_v1), the
automated *reproducible contrast* to the manual/literature ground_truth.tsv.

Dual-GT harmonization with scycle-pipeline (WORKPLAN B2, 2026-06-12): ncycle's
authoritative reference is the manual `ground_truth.tsv` (built by
build_ground_truth.py from literature + protein-level curation). This script adds
the second, fully-reproducible GT — derived purely from KEGG per-organism KO
annotation, exactly like scycle-pipeline/validation/build_ground_truth.py. For
every KO-anchored target in config/targets.yaml a genome is scored `present` iff
KEGG assigns any of that target's KOs to it (REST `link/ko/<org>`), else `absent`.
Score the pipeline against it with:  GT_FILE=kegg_ground_truth.tsv

COVERAGE (honest limitation): KEGG only annotates cultured/deposited genomes, so
this GT covers **32 of the 39 panel genomes**. The 7 excluded are precisely the
hard-to-culture N-cycle specialists with no KEGG genome: Methylacidiphilum
fumariolicum SolV, Nitrospina gracilis 3/211, comammox Nitrospira nitrosa,
anammox Brocadia sinica + Scalindua brodae, Nitrolancea hollandica Lb, and
Pseudomonas stutzeri F2a (strain not in KEGG). The manual ground_truth.tsv
covers all 39 — which is one reason it remains the authoritative reference.

CIRCULARITY (honest limitation): KEGG KO assignment is KOfam best-hit, so on the
homology-trap targets it inherits the very confound the pipeline is built to
resolve — e.g. nxrA shares K00370 with narG, so KEGG marks nxrA "present" in
denitrifiers (a GT over-call). The pipeline scored against this KEGG GT therefore
shows DEFLATED trap precision; the disagreements are the signal, reported rather
than hand-tuned away (same stance as scycle's KEGG GT). Use the manual GT for the
headline; use this one as the reproducible lower-bound contrast.

Output: validation/kegg_ground_truth.tsv (genome, target, expected, source)
Caches KEGG KO sets in validation/.kegg_cache/ (delete to refresh).
"""
from __future__ import annotations
import sys, urllib.request
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "kegg_ground_truth.tsv"
CACHE = ROOT / "validation" / ".kegg_cache"
TARGETS = ROOT / "config" / "targets.yaml"

# panel genome (proteome stem) -> KEGG organism code. All resolved via KEGG
# find/genome and validated to return >850 KO links (build session 2026-06-12).
# Exact-strain matches preferred (aha=ATCC 7966 not ML09; acp=2CP-1 not 2CP-C).
GENOME_ORG = {
    "Avinelandii_DJ": "avn",          "Neuropaea_ATCC19718": "neu",
    "Ninopinata_comammox": "nio",     # Ca. Nitrospira inopinata (comammox; IS in KEGG)
    "Njaponica_NJ11": "nja",          # KEGG "Nitrospira japonica NJ1" (same organism)
    "Nmaritimus_SCM1": "nmr",         "Nwinogradskyi_Nb255": "nwi",
    "Pdenitrificans_PD1222": "pde",   "Cnecator_H16": "reh",
    "Bsubtilis_168": "bsu",           "Synechocystis_PCC6803": "syn",
    "Ahydrophila_ATCC7966": "aha",    # exact ATCC 7966 strain
    "Kstuttgartiensis": "kst",        # Ca. Kuenenia stuttgartiensis (anammox; IS in KEGG)
    "Dvulgaris_Hildenborough": "dvu", "Scerevisiae_S288C": "sce",
    "Npcc7120": "ana",                "Smeliloti_1021": "sme",
    "Paeruginosa_PAO1": "pae",        "Wsuccinogenes_DSM1740": "wsu",
    "Ecoli_K12_MG1655": "eco",        "Soneidensis_MR1": "son",
    "Noceani_ATCC19707": "noc",       "Nmultiformis_ATCC25196": "nmu",
    "Nviennensis_EN76": "nvn",        "Hpylori_26695": "hpy",
    "Spneumoniae_ref": "spn",
    "Lacidophilus_4356": "lac",       # KEGG NCFM strain (ATCC 4356 not in KEGG); negative control
    "Aterreus_NIH2624": "ate",        # negative-control fungus
    # hold-out genomes present in KEGG (tagged holdout_v3 below so the scorer splits them)
    "Bdiazoefficiens_USDA110": "bja", "Rpalustris_CGA009": "rpa",
    "Adehalogenans_2CP1": "acp",      # exact 2CP-1 strain
    "Nhalophilus_Nc4": "nhl",         "Mcapsulatus_Bath": "mca",
}

# hold-out panel members (manual GT source=holdout_v3) that ALSO exist in KEGG.
# Tagged holdout_v3 here so score_ncycle's train/hold-out split matches the manual GT.
# (Nhollandica_Lb, Sbrodae, Sstutzeri_F2a are hold-out too but absent from KEGG.)
HOLDOUT = {"Bdiazoefficiens_USDA110", "Rpalustris_CGA009", "Adehalogenans_2CP1",
           "Nhalophilus_Nc4", "Mcapsulatus_Bath"}

# Archaeal ammonia oxidizers in the covered panel (for the KO-less amoA_archaeal
# target — PF12942-specific, no KOfam KO, so KEGG cannot assign it; derived
# explicitly, mirroring scycle's KO-less handling).
AOA_GENOMES = {"Nmaritimus_SCM1", "Nviennensis_EN76"}


def kegg_ko_set(org: str) -> set[str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{org}.ko"
    if not f.exists():
        url = f"https://rest.kegg.jp/link/ko/{org}"
        print(f"  fetching KEGG KO set: {org}", file=sys.stderr)
        with urllib.request.urlopen(url, timeout=60) as r:
            f.write_bytes(r.read())
    kos = set()
    for line in f.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) == 2 and parts[1].startswith("ko:"):
            kos.add(parts[1][3:])
    return kos


def main() -> None:
    cfg = yaml.safe_load(open(TARGETS))
    targets = cfg["targets"]
    ko_by_target = {t["id"]: list(t.get("ko") or []) for t in targets}
    ko_less = [t["id"] for t in targets if not t.get("ko")]   # amoA_archaeal

    org_kos = {g: kegg_ko_set(org) for g, org in GENOME_ORG.items()}

    rows = [("genome", "target", "expected", "source")]
    for g in GENOME_ORG:
        kos = org_kos[g]
        src = "holdout_v3" if g in HOLDOUT else "kegg_v1"
        # KO-anchored targets: present iff KEGG assigns any of the target's KOs.
        for t in targets:
            tid = t["id"]
            if tid in ko_less:
                continue
            present = any(k in kos for k in ko_by_target[tid])
            rows.append((g, tid, "present" if present else "absent", src))
        # KO-less amoA_archaeal: KEGG cannot distinguish archaeal amoA (shares
        # K10944 with bacterial amoA/pmoA), so derive explicitly from AOA membership.
        for tid in ko_less:
            present = g in AOA_GENOMES
            rows.append((g, tid, "present" if present else "absent", src))

    with open(OUT, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    n_p = sum(1 for r in rows[1:] if r[2] == "present")
    n_a = sum(1 for r in rows[1:] if r[2] == "absent")
    n_hold = len(HOLDOUT)
    print(f"[kegg_ground_truth] {len(rows)-1} cells across {len(GENOME_ORG)} genomes "
          f"({n_p} present, {n_a} absent; {n_hold} hold-out in KEGG) -> {OUT}")
    print(f"[kegg_ground_truth] NOTE: 7 panel genomes absent from KEGG (anammox / comammox / "
          f"Nitrospina / Nitrolancea / Methylacidiphilum / P. stutzeri F2a); manual GT covers all 39.")


if __name__ == "__main__":
    main()
