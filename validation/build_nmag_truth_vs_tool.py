#!/usr/bin/env python3
"""
build_nmag_truth_vs_tool.py — assemble the N-MAG realism headline table (B7),
the nitrogen analogue of scycle's p3_truth_vs_tool.tsv.

The trap: nxrA (NOB/comammox nitrite oxidation) vs narG (denitrifier nitrate
reduction) share KO K00370 + Pfam. METABOLIC and DRAM detect the shared KO/module
but cannot say which direction; ncycle resolves it (custom NOB-clade HMM + operon
synteny on nucleotide input). This table reports, per MAG: paper phenotype, the
ncycle trap call, what METABOLIC/DRAM saw (the undirected shared KO/module), and
whether the trap was resolved.

Inputs (run AFTER the comparator jobs + normalizers):
  results_nmag/<MAG>/calls/ncycle_calls.tsv        (ncycle per-target status)
  comparators/nmag/metabolic.tsv  (genome, ko)     (build_metabolic_tsv.py output)
  comparators/nmag/dram.tsv       (genome, function, present)
Output: validation/nmag_truth_vs_tool.tsv
"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results_nmag"
METABOLIC = ROOT / "comparators" / "nmag" / "metabolic.tsv"
DRAM = ROOT / "comparators" / "nmag" / "dram.tsv"
OUT = ROOT / "validation" / "nmag_truth_vs_tool.tsv"

# paper phenotype + the diagnostic gene the trap hinges on (from DOSSIER_NMAG.md)
PHENO = {
    "D1_Paracoccus_soil":          ("denitrifier",        "narG"),
    "D2_Thiobacillus_groundwater": ("denitrifier",        "narG"),
    "D3_Thioglobus_OMZ":           ("denitrifier (N2O)",  "napA"),
    "N1_Nitrosopumilus_marine":    ("AOA",                "amoA_archaeal"),
    "N2_Nitrosomonas_AOB":         ("AOB",                "amoA"),
    "N3_Nitrospina_NOB":           ("NOB",                "nxrA"),
    "N4_Nitrospira_comammox":      ("comammox",           "amoA+nxrA"),
    "A1_Brocadia_reactor":         ("anammox",            "hzsA"),
    "Neg1_Bacteroides_gut":        ("negative",           "-"),
    "Neg2_Ecoli_K12":              ("negative (assim narG)", "-"),
}
# the K00370/K00371 shared trap KOs + amo KO (what the KO-only tools collapse)
TRAP_KOS = {"K00370", "K00371"}      # nxrA/narG, nxrB/narH (undirected in KO tools)


def ncycle_calls(mag: str) -> dict[str, str]:
    f = RESULTS / mag / "calls" / "ncycle_calls.tsv"
    out = {}
    if f.exists():
        for r in csv.DictReader(open(f), delimiter="\t"):
            out[r["target_id"]] = r.get("status", "")
    return out


def load_metabolic() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    if METABOLIC.exists():
        for r in csv.DictReader(open(METABOLIC), delimiter="\t"):
            out.setdefault(r["genome"], set()).add(r["ko"])
    return out


def load_dram() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    if DRAM.exists():
        for r in csv.DictReader(open(DRAM), delimiter="\t"):
            if str(r.get("present", "")).lower() in ("1", "true", "yes"):
                out.setdefault(r["genome"], set()).add(r["function"])
    return out


def present(status: str) -> bool:
    return status in ("confirmed", "domain-only", "narrow-no-IPR")


def main() -> None:
    mtab, dram = load_metabolic(), load_dram()
    rows = [("sample", "expected_phenotype", "diagnostic_gene", "ncycle_trap_call",
             "metabolic_trap_KOs", "dram_modules", "trap_resolution")]
    for mag, (pheno, diag) in PHENO.items():
        c = ncycle_calls(mag)
        # ncycle's resolution of the nxr/nar (+amo) trap
        ncalls = [t for t in ("nxrA", "nxrB", "narG", "narH", "napA",
                              "amoA", "amoA_archaeal", "hzsA") if present(c.get(t, ""))]
        ncycle_call = "+".join(ncalls) if ncalls else "(none)"
        # what the KO-only tools saw at the shared trap loci
        mt = mtab.get(mag, set())
        mtrap = ",".join(sorted(TRAP_KOS & mt)) or "-"
        dmods = ";".join(sorted(dram.get(mag, set()))) or "-"
        # resolution verdict: did ncycle assign a direction the KO tools can't?
        has_trap_ko = bool(TRAP_KOS & mt)
        nxr = present(c.get("nxrA", "")) or present(c.get("nxrB", ""))
        nar = present(c.get("narG", "")) or present(c.get("narH", ""))
        if has_trap_ko and (nxr ^ nar):
            verdict = "✓ direction resolved (KO tools cannot)"
        elif has_trap_ko and not (nxr or nar):
            verdict = "trap KO present but ncycle non-call"
        elif pheno.startswith("negative"):
            verdict = "✓ negative" if not ncalls or ncalls == ["narG"] else "check"
        else:
            verdict = "✓" if ncalls else "miss"
        rows.append((mag, pheno, diag, ncycle_call, mtrap, dmods, verdict))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"[nmag] wrote {OUT} ({len(rows)-1} MAGs)")
    for r in rows[1:]:
        print("  " + "\t".join(r[:1]) + f"  {r[1]:<22} ncycle={r[3]:<22} {r[6]}")


if __name__ == "__main__":
    main()
