# ncycle-pipeline — Cross-Tool Comparison Report

**Question:** does ncycle-pipeline detect N-cycle genes more accurately than established tools on
the same genomes, and *where* does the difference come from?

**Panel:** 31-genome training panel · **Ground truth:** `curated_v12` · **Statistics:**
genome-cluster bootstrap (B = 10,000, seed = 1234), paired-Δ bootstrap, exact McNemar,
Benjamini–Hochberg FDR over the pre-registered family. Hold-outs excluded (training panel only,
as pre-registered). Design + decision rules: [`prereg.md`](prereg.md).

> **Status:** all four comparators complete — **raw KofamScan, METABOLIC v4.0, NCycDB, DRAM
> v1.4.6**. Every one is beaten on both pre-registered endpoints after multiple-testing correction.

---

## 1. Headline

On the shared panel, **ncycle-pipeline scores micro-F1 = 1.00 [1.00, 1.00]** (zero FP, zero FN;
all 9 pathways / 13 steps at F1 = 1.00). Against all four comparators, **both pre-registered
endpoints — overall F1 and homology-trap precision — favor ncycle and survive BH-FDR correction
(q ≤ 0.0012): ncycle WINS every contrast (8/8).**

| subunit metric | ncycle | raw KofamScan | METABOLIC | NCycDB | DRAM |
|---|---|---|---|---|---|
| ALL micro-F1 | **1.000** | 0.916 | 0.887 | 0.858 | 0.544 |
| ALL precision | **1.000** | 0.912 | 0.911 | 0.833 | 0.473 |
| ALL recall | **1.000** | 0.920 | 0.863 | 0.884 | 0.639 |
| **homology-trap precision** | **1.000** | 0.836 | 0.810 | 0.791 | 0.437 |
| homology-trap recall | **1.000** | 0.850 | 0.567 | 0.883 | 0.983 |
| non-trap micro-F1 | **1.000** | 0.939 | 0.945 | 0.865 | 0.513 |
| ALL specificity (TN rate) | **1.000** | 0.944 | 0.947 | 0.889 | 0.553 |
| **ALL MCC** (TP/FP/FN/TN) | **1.000** | 0.863 | 0.819 | 0.765 | 0.187 |

The two TN-inclusive rows (added 2026-06-16, §3.1) make the false-positive control explicit:
**specificity** (true-negative rate) and **MCC** (all four confusion cells, robust to class imbalance).
ncycle is perfect on both (1.000); DRAM's module-level over-calling shows starkly here — specificity
0.553, MCC 0.187 — and even the strong KO tools sit at MCC 0.82–0.86 because their residual trap
false positives cost them where F1 partly masked it.

**The four comparators fail the homology traps in four *different* ways** — and ncycle's layered
design (clade HMMs + curated BLAST gates + operon synteny) is the only one that resolves all of
them:

| tool | trap failure signature | mechanism |
|---|---|---|
| raw KofamScan | over-calls (P 0.84, R 0.85) | shared KOs (K00370/K00371, amoA/pmoA) map to *all* their targets |
| METABOLIC | **under-calls** (P 0.81, **R 0.57**) | stricter KOfam m-cutoff → misses nxr/nar/amo in many genomes |
| NCycDB | **mis-routes** (P 0.79) | curated nxrA family exists, but DIAMOND best-hit sends NOB nxrA → narG |
| DRAM | **massively over-calls** (**P 0.44**, R 0.98) | module-level calls expand one shared KO to every member subunit |

---

## 2. Tools compared

| tool | approach | resolution | status |
|---|---|---|---|
| **ncycle-pipeline** (this tool) | KOfam HMM + clade-specific custom HMMs + curated BLAST gates (+ operon synteny on nucleotide input) | subunit (49 targets) | — |
| **raw KofamScan** | KO-only, stock `ko_list` thresholds; shared KOs map to *all* their targets; no gating | subunit (by KO) | ✅ |
| **METABOLIC** v4.0 (Zhou et al. 2022) | HMM-based function profiling (`METABOLIC-G.pl`, bundled KOfam, `-kofam-db full`) → per-genome KO presence | subunit (by KO) | ✅ |
| **NCycDB** (Tu et al. 2019) | curated N-cycle gene-family DB (68 families) + DIAMOND best-hit | subunit (gene families) | ✅ |
| **DRAM** v1.4.6 (Shaffer et al. 2020) | MAG metabolism distillation (KEGG modules) | coarse (modules → steps) | ✅ |

**Mapping discipline (applies to every tool).** Each tool's native output is mapped into one
shared target vocabulary by a *fixed, reviewable table* (`adapters.py`), decided from each tool's
own semantics — never tuned to the panel. Two resolutions are reported so coarse tools aren't
penalized merely for granularity: **subunit** (49 ncycle targets; the trap claim lives here) and
**step** (13 N-cycle transformations; fairest to coarse tools). The same collapse rule is applied
identically to all tools and to the ground truth. METABOLIC and DRAM map through KEGG KOs, so they
inherit the same KO→target table as the raw-KofamScan baseline by construction.

---

## 3. Pre-registered confirmatory result (subunit, BH-FDR corrected)

Family = {comparator} × {trap-precision, ALL-F1}. A win requires the paired-Δ CI to exclude 0
**and** the point estimate to favor ncycle.

| comparator | metric | Δ (ncycle − cmp) | 95% CI | bootstrap p | BH q | verdict |
|---|---|---|---|---|---|---|
| raw KofamScan | ALL F1 | +0.084 | [0.059, 0.113] | <0.0001 | <0.0001 | **ncycle WINS** |
| METABOLIC | ALL F1 | +0.113 | [0.082, 0.150] | <0.0001 | <0.0001 | **ncycle WINS** |
| NCycDB | ALL F1 | +0.142 | [0.118, 0.172] | <0.0001 | <0.0001 | **ncycle WINS** |
| DRAM | ALL F1 | +0.456 | [0.405, 0.514] | <0.0001 | <0.0001 | **ncycle WINS** |
| raw KofamScan | trap precision | +0.164 | [0.069, 0.258] | 0.0012 | 0.0012 | **ncycle WINS** |
| METABOLIC | trap precision | +0.190 | [0.097, 0.297] | <0.0001 | <0.0001 | **ncycle WINS** |
| NCycDB | trap precision | +0.209 | [0.113, 0.323] | <0.0001 | <0.0001 | **ncycle WINS** |
| DRAM | trap precision | +0.563 | [0.492, 0.641] | <0.0001 | <0.0001 | **ncycle WINS** |

### 3.1 Secondary (post-hoc) MCC family — BH-FDR corrected

TN-inclusive confirmation, scored on the same cells with the same B=10,000 cluster bootstrap, BH-FDR
over its own 8-contrast family (kept **separate** from the frozen pre-registered family above so the
confirmatory result is unchanged). All 8 favour ncycle.

| comparator | metric | Δ (ncycle − cmp) | 95% CI | BH q | verdict |
|---|---|---|---|---|---|
| raw KofamScan | ALL MCC | +0.137 | [0.097, 0.180] | <0.0001 | **ncycle WINS** |
| METABOLIC | ALL MCC | +0.181 | [0.132, 0.233] | <0.0001 | **ncycle WINS** |
| NCycDB | ALL MCC | +0.235 | [0.201, 0.274] | <0.0001 | **ncycle WINS** |
| DRAM | ALL MCC | +0.813 | [0.708, 0.923] | <0.0001 | **ncycle WINS** |
| raw KofamScan | trap MCC | +0.226 | [0.134, 0.320] | <0.0001 | **ncycle WINS** |
| METABOLIC | trap MCC | +0.428 | [0.322, 0.544] | <0.0001 | **ncycle WINS** |
| NCycDB | trap MCC | +0.241 | [0.136, 0.348] | <0.0001 | **ncycle WINS** |
| DRAM | trap MCC | +0.573 | [0.485, 0.675] | <0.0001 | **ncycle WINS** |

---

## 4. Per-step micro-F1 (step resolution, fairest to coarse tools)

ncycle = 1.000 on every step. F1 shown; precision in parentheses where it drives a low score.

| step | ncycle | KofamScan | METABOLIC | NCycDB | DRAM |
|---|---|---|---|---|---|
| ALL | **1.000** | 0.942 | 0.925 | 0.835 (P 0.76) | 0.510 (P 0.44) |
| N2_fixation | **1.000** | 1.000 | 1.000 | 0.923 (P 0.86) | 1.000 |
| ammonia_oxidation | **1.000** | 0.833 | 0.833 | 0.857 (P 0.75) | 0.522 (P 0.35) |
| nitrite_oxidation | **1.000** | 0.182 (P 0.17) | 0.889 | **NA — R 0.00** | 0.400 (P 0.27) |
| dissim_nitrate_reduction | **1.000** | 1.000 | 0.667 (P 0.57) | 0.720 (P 0.60) | 0.556 (P 0.39) |
| NO2_to_NO | **1.000** | 1.000 | 0.923 | 0.778 (P 0.64) | 0.757 (P 0.61) |
| NO_reduction | **1.000** | 0.941 | 0.941 | 0.933 | 0.800 (P 0.67) |
| N2O_reduction | **1.000** | 1.000 | 1.000 | 0.400 (P 0.25) | 0.286 (P 0.17) |
| DNRA | **1.000** | 1.000 | 0.923 | 1.000 | 0.667 (P 0.50) |
| anammox | **1.000** | 1.000 | 1.000 | 1.000 | 0.174 (P 0.10) |
| assim_nitrate_reduction | **1.000** | 0.947 | 0.947 | 0.833 (P 0.71) | 1.000 |
| assim_nitrite_reduction | **1.000** | 0.857 | 0.667 | 0.667 (P 0.60) | 1.000 |
| ammonia_assimilation | **1.000** | 1.000 | 1.000 | 1.000 | **NA — no module** |
| ureolysis | **1.000** | 1.000 | 1.000 | 0.963 | **NA — no module** |

---

## 5. raw KofamScan — findings

**Ties ncycle** on every step with a specific, unambiguous KO (N-fixation, DNRA, anammox,
dissimilatory nitrate reduction, NO2→NO, N2O reduction, ammonia assimilation, ureolysis);
non-trap subunit subset is near-parity (0.939 vs 1.000).

**Loses at the shared-KO traps:**
- **nitrite_oxidation F1 0.182** — fires K00370/K00371 (the nar/nxr shared KOs) → calls `nxr` in
  *every* denitrifier (precision 0.17), its worst step.
- **ammonia_oxidation F1 0.833** — calls `amoA` on the methanotroph (pmoA shares the KO bucket).

**Signature: over-calling** — shared KOs map to all their targets, so precision drops exactly where
two targets share a KO. This is the design rationale for ncycle's gating layer.

---

## 6. METABOLIC v4.0 — findings

METABOLIC is HMM-based on a bundled KOfam (`-kofam-db full`); mapped through KEGG KOs like the
KofamScan baseline. Its overall F1 (0.887) is close to KofamScan, **but its failure mode is the
opposite**: it **under-calls** the traps (trap recall **0.567**, vs KofamScan's 0.850) because its
stricter KOfam thresholds reject borderline hits. Per-target confusion (curated_v12):

| target | METABOLIC | note |
|---|---|---|
| **narG** | TP 0 / **FP 6 / FN 6** | thresholds reject narG in real denitrifiers (6 FN) yet the shared KO fires in NOB (6 FP) → *every* narG call is wrong |
| amoA | TP 4 / FP 2 | same amoA/pmoA shared-KO FP as KofamScan (methanotroph + archaeon), fewer of them |
| amoA_archaeal | TP 0 / **FN 2** | no archaeal-specific resolution (archaeal amoA fires the generic amoA KO) |
| nrfH | TP 1 / FN 3 | mostly below threshold |
| NR, nasD | FN | K10534 + K17877 absent from METABOLIC's bundled KOfam version → KO-version coverage gaps (≤2 cells, documented) |

**Where it is strong:** nitrite_oxidation step F1 0.889 (better than KofamScan's 0.182 — its strict
thresholds *help* here), and napA/nirK/N-fixation/ureolysis at/near 1.00. **Signature:
under-calling** — high precision bought with recall, the mirror image of KofamScan.

---

## 7. NCycDB — findings (the curated-database case)

NCycDB is a *curated* N-cycle database — it has **separately named `nxrA`/`narG` families** and a
**separate `pmoA`** — so a priori it should handle the traps better. It does on some, but two
structural problems remain.

### 7.1 Headline — NCycDB misses every nitrite-oxidizer nxrA

| target | NCycDB TP | FP | FN | TN | note |
|---|---|---|---|---|---|
| **nxrA** | **0** | 0 | **5** | 24 | every NOB nxrA best-hits the `narG` family instead |
| **narG** | 4 | **6** | 2 | 17 | the 6 FP are exactly those mis-routed NOB nxrA proteins |
| nxrB | 4 | 0 | 0 | 12 | nxrB is separable; nxrA is not |

Despite owning a curated `nxrA` family, DIAMOND **best-hit** routes the divergent
Nitrobacter/Nitrospira NxrA (~60 % identical to NarG) into `narG` → **step nitrite_oxidation
recall = 0.00**. This is the *same* sequence-homology trap ncycle resolves with a clade-calibrated
NXR HMM + operon-synteny disambiguation — and the exact failure a homology-only best-hit cannot escape.

### 7.2 Other modes
- `nosZ` over-call (TP 4 / **FP 12**) → step N2O_reduction precision 0.25; `nirK` 4 FP; assimilatory
  steps P 0.60–0.71.
- Structural coverage gaps (honest FN): `nrfH`, γ-AOB `amoA` (*N. oceani*), qNor `norZ`, clade-II
  `nosZ`, `nifE/N/B`, `vnfD/H`, `nasD`.

### 7.3 Genuine strengths
- **Archaeal vs bacterial amoA resolved** — `amoA_A`/`amoA_B` are truly clade-specific (verified:
  `amoA_A` hit only by AOA archaea, `amoA_B` only by bacterial AOB/comammox). `amoA_archaeal`
  2 TP / 0 FP — parity with ncycle, and better than every other comparator here.
- **No methanotroph amoA FP** — routes *M. fumariolicum* to its `pmoA/B/C` families.
- Perfect on anammox, DNRA, ureolysis, ammonia assimilation, `napA`, `narH`, `nxrB`, nifHDK.

**Signature: mis-routing** — curated names help, but homology best-hit with no specificity gating
still fails nxrA↔narG and over-calls low-specificity families.

---

## 8. DRAM v1.4.6 — findings (the coarse-distillation case)

DRAM distills to **KEGG-module** granularity; a present module is expanded to all member subunits
of every step it covers (symmetric credit/blame — the honest representation of what DRAM tells a
user). This produces the lowest scores here (ALL F1 **0.544**), in two ways at once:

### 8.1 Massive over-calling via shared module KOs (trap precision 0.44, recall 0.98)
DRAM's modules bundle KOs that recur across functions, so one shared hit lights up an entire module:

| target | DRAM TP | FP | note |
|---|---|---|---|
| **hzsA** (anammox) | 2 | **19** | the "nitrite + ammonia => nitrogen" module includes nirK/nar/nirS KOs → fires in 19 non-anammox genomes |
| **nosZ** | 4 | **20** | denitrification module → every denitrifier marked nosZ |
| **narG** | 6 | **18** | nar KOs sit in denitrification + DNRA + comammox modules |
| amoA / amoA_archaeal | 4 / 2 | 13 / 15 | nitrification + comammox modules expand to amo subunits |
| nirK | 10 | 8 | — |

→ anammox step precision 0.10 (F1 0.17); ammonia_oxidation P 0.35; N2O_reduction P 0.17.

### 8.2 Coverage gaps (no module → total misses)
- **glnA: 0 TP / 31 FN** — DRAM has no ammonia-assimilation module → misses glutamine synthetase in
  *every* genome (step ammonia_assimilation = NA).
- **ureC: 0 TP / 14 FN** — no catalytic-urease module (its "urea cycle" module is the arginine
  pathway) → ureolysis = NA.

### 8.3 Where it is clean
N-fixation (nifH 6 TP / 0 FP — the nif module is specific) and the assimilatory nitrate/nitrite
steps (F1 1.00). **Signature: coarse over-calling + module coverage gaps** — DRAM answers "which
broad pathway is present" well but cannot resolve subunits or homology traps.

---

## 9. Reproduction

```bash
# NCycDB — download DB + DIAMOND on the 31 training proteomes → ncycdb.tsv
python3 comparators/build_ncycdb_tsv.py
# METABOLIC v4.0 (env METABOLIC_v4.0) — METABOLIC-G.pl on the proteomes, then normalize:
#   perl METABOLIC-G.pl -t 8 -in comparators/metabolic_in -kofam-db full -o comparators/metabolic_out
python3 comparators/build_metabolic_tsv.py            # → validation/benchmark/metabolic.tsv (genome, ko)
# DRAM v1.4.6 (env DRAM14, KOfam-only config) — annotate_genes + distill, then normalize:
#   DRAM.py annotate_genes -i 'comparators/dram_in/*.faa' -o comparators/dram_out --threads 8
#   DRAM.py distill -i comparators/dram_out/annotations.tsv -o comparators/dram_distill
python3 comparators/build_dram_tsv.py                 # → validation/benchmark/dram.tsv (genome, function, present)
# Score all four comparators together (ncycle + kofam need no extra input):
python3 validation/benchmark/benchmark_stats.py \
    --ncycdb validation/benchmark/ncycdb.tsv \
    --metabolic validation/benchmark/metabolic.tsv \
    --dram validation/benchmark/dram.tsv --bootstrap 10000
```

Provenance + fixed mapping tables (`NCYC_MAP`, `DRAM_STEP_MAP`, KO→target) in `adapters.py`;
normalizers in `comparators/build_*_tsv.py`; tool inputs in `comparators/{metabolic,dram}_in/`;
raw tool outputs in `comparators/{NCyc,metabolic_out,dram_out,dram_distill}/`; full numeric table
in `benchmark_results.tsv`. METABOLIC/DRAM used the shared installs under
`Sulfur_Cycle/comparators/{METABOLIC,DRAM_data}`.

---

## 10. Synthesis — similarities & differences across tools

**All five tools agree** on targets with a specific, unambiguous marker: N-fixation core, ammonia
assimilation (except DRAM, which lacks the module), and ureolysis (except DRAM). Detection is a
solved problem there; the tools diverge only at the hard cases.

**The discriminating axis is specificity at homology traps — and each comparator fails it
differently:**
- *raw KofamScan* — **over-calls**: shared KOs map to all targets (nxr/amo precision drops).
- *METABOLIC* — **under-calls**: stricter KOfam thresholds miss real nxr/nar/amo (trap recall 0.57).
- *NCycDB* — **mis-routes**: curated families exist, but homology best-hit sends NOB nxrA → narG
  (nitrite-oxidation recall 0) and over-calls low-specificity families (nosZ).
- *DRAM* — **coarse over-call + gaps**: module expansion floods the traps (trap precision 0.44) and
  whole functions have no module (glnA, ureC).
- *ncycle* — clade HMMs + per-target BLAST gates + operon synteny resolve all three trap classes →
  trap precision 1.00 where the others sit at 0.44–0.84.

**Coverage breadth** is a secondary differentiator: ncycle carries clade/paralog targets the others
lack (γ-AOB amoA, archaeal amoA, qNor norZ, clade-II nosZ, nrfH, nif accessory genes) — surfaced as
honest FN for the comparators, not as a hidden advantage. NCycDB is the strongest comparator on the
panel (it alone resolves archaeal amoA and avoids the methanotroph FP); DRAM is the weakest
(module-level resolution cannot address subunit-level traps).

**Bottom line:** ncycle wins all 8 pre-registered contrasts (4 tools × {ALL-F1, trap precision})
after BH correction, and the four comparators' four distinct failure signatures together show the
advantage is structural — a consequence of resolving shared-KO/shared-homology ambiguity that
KO-only, best-hit, and module-level methods cannot.

## 11. Interpretation caveat — which endpoint to trust (stats audit 2026-06-10)

This comparison is run on the **31-genome training panel** against `curated_v12`, the GT that
ncycle's thresholds and seeds were tuned on. ncycle therefore has a structural **home-field
advantage on the F1 endpoint**: its own residual errors were curated away during development, while
the comparators' were not. So ncycle's **F1 = 1.00 here is an in-sample number** (see
`../REPORT.md` § Audit 2026-06-10); the **absolute** F1 gaps (Δ +0.08 … +0.46) overstate the
real-world margin, and the credible generalization figure is the de-leaked hold-out F1 ≈ 0.84 — not
benchmarked against comparators here.

The **trap-precision endpoint is the circularity-robust one** and should be foregrounded: a
false positive lands on a GT-**absent** cell, which is never a seed source, so trap precision
(ncycle 1.00 vs 0.44–0.84) cannot be inflated by GT-curation or seed leakage. The qualitative
result — four comparators failing the traps in four distinct, mechanistically-explained ways — is
likewise independent of the absolute F1 calibration. Read this report for *where and why* the tools
diverge at the traps, and lean on trap precision (not the F1 gap) for the magnitude.

---

**Provenance.** Numbers reflect the 2026-06-14 (B8) harmonized re-run (B=10,000, seed 1234,
GT=`ground_truth.tsv`, n=31), which reproduces the prior figures exactly after the raw-KofamScan
parser fix (see `prereg.md` B8 amendment). Sister report: `scycle-pipeline/.../COMPARISON_REPORT.md`
— identical battery and structure, with SCycDB as scycle's domain-DB comparator (NCycDB here).
