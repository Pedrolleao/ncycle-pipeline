# ncycle-pipeline — Validation Report

**Tool:** ncycle-pipeline (nitrogen-cycle pathway/protein detection from genomes, MAGs, and proteomes)
**Date:** 2026-05-31 · **Ground truth:** `curated_v12` (+ `holdout_v3`, `expanded_tn_v12`) · **Status:** publication-ready (Phase 5 complete)

This report documents *how* the method was validated and the *measured values* of everything
tested. It is the consolidated summary of the running log in `validation/REPORT.md` (see
§ Status + § Audit 2026-05-30 there for full provenance); reproduction commands are at the
end. Headline: **training micro-F1 1.00 [1.00, 1.00]** (zero FP, zero FN; all 9 pathways at
F1 = 1.00) on a 31-genome training panel, with an **honest 8-genome / 5-phyla hold-out
F1 = 0.85 [0.70, 0.95]** as the credible generalization result, and homology-trap precision
1.00 — beating raw KofamScan (F1 1.00 vs 0.93; FP 0 vs 14).

> **Provenance note.** An earlier version of this file reported the curated_v4 era
> (training F1 0.97, 2-genome hold-out F1 1.00). Those numbers were superseded by the
> 2026-05-30 dual-agent audit and the Phase 5 work that followed: a publication-blocking
> panel contamination was fixed (P5.0), bootstrap CIs were wired into the headline path
> (P5.1), and the hold-out was expanded from 2 to 8 genomes spanning 5 phyla (P5.3.5),
> which gave the hold-out CI real empirical width. The training F1 rose to 1.00 (verified,
> zero residual FP/FN); the hold-out F1 settled at 0.85 because the expanded hold-out
> exposes genuine cross-phylum classification heterogeneity that the old 2-genome hold-out
> could not. The 0.85 is the number to cite for generalization.

---

## 1. Headline results

| Metric | Value |
|---|---|
| **Training micro-F1** (31 genomes, 645 cells) | **1.00 [1.00, 1.00]** |
| Training micro-precision / recall | 1.00 / 1.00 |
| Training TP / FP / FN / TN | 249 / **0** / **0** / 396 |
| Pathways at F1 = 1.00 | **9 / 9** |
| Homology-trap precision | **1.00 [1.00, 1.00]** |
| **Hold-out micro-F1** (8 genomes, 5 phyla, 194 cells) | **0.85 [0.70, 0.95]** |
| Hold-out precision / recall | 0.86 / 0.83 |
| Hold-out TP / FP / FN / TN | 75 / 12 / 15 / 92 |
| vs raw KofamScan (same panel) | **F1 1.00 vs 0.93; FP 0 vs 14** |

All four "robust" acceptance criteria (coverage, accuracy-with-CIs, traps-separated,
benchmark) are met — see § 6.

---

## 2. Validation methodology

### 2.1 Reference panel (39 genomes: 31 training + 8 hold-out)
A curated panel spanning all 9 N-cycle pathways with ≥3 positive genomes per major pathway,
homology-trap decoys, and clean negatives. Split into a **training panel (31 genomes)** used
for scoring/tuning and an **8-genome independent hold-out** drawn from genera/phyla absent
from every HMM/BLAST training set.

- **Training (31):** Ahydrophila_ATCC7966, Aterreus_NIH2624, Avinelandii_DJ, Bsinica_JPN1,
  Bsubtilis_168, Cnecator_H16, Dvulgaris_Hildenborough, Ecoli_K12_MG1655, Hpylori_26695,
  Kstuttgartiensis, Lacidophilus_4356, Mfumariolicum_SolV, Neuropaea_ATCC19718, Ngracilis_3211,
  Ninopinata_comammox, Njaponica_NJ11, Nmaritimus_SCM1, Nmultiformis_ATCC25196, Nnitrosa_comammox,
  Noceani_ATCC19707, Npcc7120, Nviennensis_EN76, Nwinogradskyi_Nb255, Paeruginosa_PAO1,
  Pdenitrificans_PD1222, Scerevisiae_S288C, Smeliloti_1021, Soneidensis_MR1, Spneumoniae_ref,
  Synechocystis_PCC6803, Wsuccinogenes_DSM1740.
- **Hold-out (8), spanning 5 phyla (α-/δ-/γ-proteo, Chloroflexi, Planctomycetes):**
  Bdiazoefficiens_USDA110 (*Bradyrhizobium*), Rpalustris_CGA009 (*Rhodopseudomonas*),
  Adehalogenans_2CP1 (*Anaeromyxobacter*, δ-proteo, clade-II NosZ carrier),
  Nhalophilus_Nc4 (*Nitrosococcus*, γ-AOB), Nhollandica_Lb (*Nitrolancea*, Chloroflexi NOB),
  Sbrodae (*Ca. Scalindua*, Planctomycetes anammox), Sstutzeri_F2a (γ-proteo denitrifier),
  Mcapsulatus_Bath (*Methylococcus*, γ-proteo methanotroph — amoA decoy).

Decoys deliberately included to test the homology traps: *Methylacidiphilum* / *Methylococcus*
(methanotroph; pmoA vs amoA), *Desulfovibrio* (sulfate reducer; dsrAB vs nir), *S. cerevisiae*
(eukaryote negative), and both NOB and denitrifiers (NxrA vs NarG).

### 2.2 Panel QC gate (`validation/validate_panel.py`)
After P5.0 caught a contaminated panel file (a "*Nitrobacter*" proteome that was actually
*Burkholderia thailandensis*, zero Nitrobacter tags), a per-genome QC gate was added: it
samples protein headers and fails fast if <60% match the expected genus/species for the
filename. Wired into the Makefile (`validate-panel`) as a regression dependency, so a
wrong-genus swap can never silently reach scoring again.

### 2.3 Ground truth (`curated_v12`)
Literature/annotation-based present/absent calls, KEGG per-organism KO-verified, with an
environmental-microbial-genomics specialist audit. **Every genome's identity was verified by
gene content, not by label** (this discipline caught multiple accession/identity bugs during
construction). GT was regenerated deterministically and versioned through `curated_v12`;
v11 resolved the P5.0 contamination (narG stands, nosZ + norZ corrections retracted), and
v12 added a programmatic true-negative expansion for the two new clade HMMs (amoA_gamma,
nosZ_clade2).

- **839 total cells** across 39 genomes: source tags `curated_v11` (training, manually
  curated), `holdout_v3` (8 hold-out genomes), `expanded_tn_v12` (programmatic TN expansion).
- Training scoring = **645 cells** (31 genomes); hold-out scoring = **194 cells** (8 genomes).
- The scorer scores **only curated cells**; unlisted (genome, target) pairs are not scored.

### 2.4 Scoring
Cell-level binary classification (present vs absent); positive class = "present". Matrix status
codes map: 2 confirmed / 1 domain-only or narrow-no-IPR → **present**; 0 absent / −1 disqualified
→ **absent**. Metrics: per-target and micro-aggregated precision, recall, F1; median per-target
F1; reported overall, on the **homology-trap subset**, the **non-trap subset**, and per-pathway.
Every aggregate, per-pathway, and hold-out estimate is reported with a 95% genome-cluster
bootstrap CI (§ 4).

### 2.5 Validation strategies applied
1. **Held-out cross-target scoring** against curated ground truth (the core accuracy measurement).
2. **Independent hold-out genomes** (8 genomes / 5 phyla, unseen genera) — overfitting test.
3. **Leave-one-genus-out cross-validation** of the custom HMMs — generalization test.
4. **Comparator benchmark** vs raw KofamScan, with genome-level bootstrap CIs, paired
   significance tests, and multiple-testing correction (harness in `validation/benchmark/`).
5. **Independent-modality validation** of the operon-synteny trap resolver on nucleotide input.
6. **Regression gate** with CI-aware floors locking the metrics against silent degradation.
7. **Unit tests** of the statistical functions (McNemar, BH, paired bootstrap).
8. **Provenance pinning** — KOfam release sha256-pinned, UniProt seed snapshot cache +
   Pfam pin, `verify_seeds.py` drift check.

---

## 3. Results — all measured values

### 3.1 Aggregate metrics

```
TRAINING (31 genomes, 645 cells):
  TP=249  FP=0  FN=0  TN=396
  micro-P = 1.00 [1.00, 1.00]   micro-R = 1.00 [1.00, 1.00]   micro-F1 = 1.00 [1.00, 1.00]
  homology-trap precision = 1.00 [1.00, 1.00]

HOLD-OUT (8 genomes, 5 phyla, 194 cells):
  TP=75  FP=12  FN=15  TN=92
  P = 0.86 [0.75, 0.95]   R = 0.83 [0.64, 0.98]   F1 = 0.85 [0.70, 0.95]
```

The training CIs collapse to [1.00, 1.00] because the panel has zero FP and zero FN globally
(honest structural behavior; a synthetic-imperfect-panel sanity check confirms the bootstrap
function produces real width on imperfect data). The hold-out CI has genuine empirical width.

### 3.2 Per-pathway micro-F1 (training panel)

All 9 pathways score F1 = 1.00 [1.00, 1.00] (645-cell panel, FP=0/FN=0):
nitrogen_fixation, ammonia_assimilation, assimilatory, anammox, nitrification_ammonia,
nitrification_nitrite, denitrification, dnra, organic_n_mineralization.

### 3.3 Per-target metrics (training panel)
Of the 49 configured targets carrying curated cells, **every scored target is at
precision = recall = F1 = 1.00** on the training panel (full table in
`validation/ncycle_metrics.tsv`). `nosZ_clade2` has no training positive (n_pos = 0 → F1 NA;
its first real positive is in the hold-out, *Anaeromyxobacter dehalogenans*). Per-target rows
carry `n_pos` (TP+FN) and an `f1_reliable` flag (true when n_pos ≥ 3); low-positive targets
are flagged honestly rather than counted as strong evidence.

### 3.4 Independent hold-out test (generalization)
The 8-genome / 5-phyla hold-out is the headline generalization result: **F1 = 0.85
[0.70, 0.95]**. The 12 FP + 15 FN are documented as honest findings and **deliberately NOT
back-edited into the GT** (which would convert the hold-out into a training set). They are
inventoried in § 7. amoA_gamma and nosZ_clade2 gain their first hold-out positives here.

### 3.5 Leave-one-genus-out cross-validation (custom HMMs)
- **nxrA** trained without the *Nitrospira* genus (only *Nitrobacter* + *Nitrospina*) still
  scores held-out *Nitrospira* genomes at **985–996 ≫ TC 700** → generalizes across genus.
- **amoA** trained AOB-only scores comammox at **≈TC 400 (borderline)** → amoA is
  clade-sensitive; including comammox seeds is a real coverage need, not overfitting.

### 3.6 Comparator benchmark (vs four external tools)
Same training panel (31 genomes), same ground truth (`curated_v12`), same scorer; genome-cluster
bootstrap (B=10,000), paired-Δ, exact McNemar, BH-FDR over the pre-registered family
{comparator}×{trap-precision, ALL-F1}. Pre-registration in `validation/benchmark/prereg.md`; full
narrative in `validation/benchmark/COMPARISON_REPORT.md`.

| subunit metric | ncycle | KofamScan | METABOLIC | NCycDB | DRAM |
|---|---|---|---|---|---|
| ALL micro-F1 | **1.000** | 0.916 | 0.887 | 0.858 | 0.544 |
| homology-trap precision | **1.000** | 0.836 | 0.810 | 0.791 | 0.437 |
| homology-trap recall | **1.000** | 0.850 | 0.567 | 0.883 | 0.983 |

**ncycle beats all four comparators on both pre-registered endpoints — 8/8 contrasts win after
BH-FDR (q ≤ 0.0012).** Each comparator fails the homology traps in a different way, which is the
core finding:
- **raw KofamScan — over-calls:** shared KOs (K00370/K00371; amoA/pmoA) map to *all* their targets
  → nxr called in every denitrifier (nitrite-ox F1 0.18), amoA on the methanotroph.
- **METABOLIC v4.0 — under-calls:** stricter KOfam thresholds reject borderline hits → trap recall
  0.567; e.g. narG = 0 TP / 6 FP / 6 FN. (KO-version gaps: lacks K10534/K17877 = NR/nasD.)
- **NCycDB — mis-routes:** despite a curated `nxrA` family, DIAMOND best-hit sends every NOB nxrA
  into `narG` (nxrA 0 TP / 5 FN; nitrite-ox recall 0; 6 narG FP); over-calls nosZ (12 FP). *But*
  it alone resolves archaeal vs bacterial amoA (2 TP / 0 FP) and avoids the methanotroph FP.
- **DRAM v1.4.6 — coarse over-call + gaps:** module-level calls expand one shared KO to all member
  subunits (anammox hzsA 2 TP / 19 FP; nosZ 4 / 20; trap precision 0.44), and it has no module for
  glnA (0/31) or catalytic ureC (0/14).

ncycle's clade HMMs + curated BLAST gates + operon synteny resolve all three trap classes → trap
precision 1.00 where the comparators sit at 0.44–0.84. Provenance + fixed mapping tables in
`adapters.py`; runners `comparators/build_{ncycdb,metabolic,dram}_tsv.py`.

### 3.7 Operon-synteny trap resolution — independent nucleotide-mode validation
The one sequence-inseparable case (*Nitrobacter* NxrA ≈ NarG, ~60% identity, same score band)
is resolved by gene-neighborhood context on nucleotide/MAG input (narGHJI operon synteny vs
standalone NXR locus). Validated end-to-end (Prodigal → synteny) on *N. winogradskyi* Nb-255
(recovers the proteome-mode nxrA/nxrB call) and *E. coli* K-12 (operon-only → narG, no nxr
false-positive). On pre-called proteomes without gene coordinates this case remains a
documented limit, not a bug.

### 3.8 Validation progression across hardening phases

| Phase | Panel / GT | training F1 | trap | FP | note |
|---|---|---|---|---|---|
| P1 baseline | 15 / curated_v1 | 0.84 | 0.54 | 14 | KO + auto-fetched seeds; deficit isolated to traps |
| P2 curated seeds | 15 / curated_v1 | 0.94 | 0.80 | 2 | curated clade BLAST seeds + diamond-DB bugfix |
| P3 custom HMMs | 16 / curated_v2 | 0.95 | 0.88 | 3 | nxrA/nxrB/amoA clade HMMs; comammox validated |
| Panel expansion | 33 / curated_v3 | 0.95 | 0.86 | — | +17 genomes, +2 hold-outs |
| P4 fixes + CV | 33 / curated_v4 | 0.97 | 0.96 | 2 | dropped harmful gates; beats KofamScan; LOGO-CV |
| **P5 audit response** | **31 + 8 / curated_v12** | **1.00 [1.00, 1.00]** | **1.00** | **0** | contamination fix, bootstrap CIs, 8-genome hold-out (F1 0.85) |

### 3.9 Regression gate (CI-aware floors)
`validation/test_regression.py` — fails (exit 1) if any metric's CI lower bound drops below a
pinned v12 floor. Current run: **all checks PASS.**

| Metric | v12 observed CI | Floor (CI lower bound) |
|---|---|---|
| ALL micro-F1 | [1.00, 1.00] | ≥ 1.00 |
| ALL precision | [1.00, 1.00] | ≥ 1.00 |
| trap precision | [1.00, 1.00] | ≥ 0.95 |
| hold-out F1 | [0.695, 0.948] | ≥ 0.65 |
| per-pathway F1 | all [1.00, 1.00] | ≥ 0.70 |
| ALL false positives (point check) | 0 | ≤ 4 |
| scored cells (point check) | 645 | ≥ 600 |

---

## 4. Statistical methods
- **Unit of independence = genome** (cells within a genome are correlated); all 95% CIs use a
  **genome-level cluster bootstrap** (resample genomes with replacement, B = 10,000, seed = 1234,
  percentile method per `validation/benchmark/benchmark_stats.py:130-144`), avoiding pseudoreplication.
- **Paired significance:** bootstrap of Δ = metric(ncycle) − metric(comparator) on the *same*
  resampled genomes; a win requires the Δ CI to exclude 0 (two-sided bootstrap p reported).
- **Secondary:** exact (binomial) McNemar on cell-level discordant pairs.
- **Multiple testing:** Benjamini–Hochberg FDR over the pre-specified primary family.
- **GT-vs-detection attribution:** `validation/decompose_gt.py` walks documented GT cell-flips
  in reverse and re-scores; holding detection constant, GT revisions account for only +0.016 F1
  — most of the gain from the P1 baseline (F1 ≈ 0.84) is detection-side.
- Statistical functions unit-tested against known values (McNemar b=10,c=0 → p=0.00195;
  symmetric → p=1; BH and paired self-Δ checks pass).

---

## 5. Reproducibility
- **KOfam release pinned by sha256** (`b03d20b9…`, 2026-05-24); the build verifies the tarball
  and writes a provenance stamp (`resources/.cache/kofam_release.txt`).
- **UniProt seed snapshot cache + Pfam-A sha256 pin** (124 accessions across 32 targets +
  PF12942); `validation/verify_seeds.py` (`make verify-seeds`) catches silent UniProt revisions
  before they can invalidate TC calibration.
- Standalone conda env (`envs/ncycle.yaml`); `run.py` defaults to the `ncycle-pipeline` env
  (`NCYCLE_ENV` to override).
- Ground truth regenerated deterministically from `validation/build_ground_truth.py`.

```bash
python3 validation/build_ground_truth.py                 # writes ground_truth.tsv (curated_v12)
python3 validation/validate_panel.py                     # panel-QC gate (genus/species check)
python run.py --input ../test_panel --skip-db-setup      # proteome panel → results/
python3 validation/score_ncycle.py                       # §3 tables, with bootstrap CIs
python3 validation/compare_kofam.py                      # §3.6 raw-KofamScan comparator
python3 validation/test_regression.py                    # §3.9 CI-aware gate
```

---

## 6. Acceptance criteria — all met

1. **Coverage** — every pathway has ≥1 training TP; every clade has ≥1 hold-out TP after the
   P5.3.5 expansion (incl. first real positives for amoA_gamma and nosZ_clade2). ✅
2. **Accuracy** — F1 reported with 95% genome-cluster bootstrap CIs; median F1 ≥ 0.90
   (training 1.00; hold-out 0.85 ± ~0.13). ✅
3. **Traps separated** — amoA↔pmoA, nxr↔nar, nirB/D/A↔dsr all clade-discriminated; training
   trap precision = 1.00 [1.00, 1.00]. ✅
4. **Benchmark** — pipeline vs raw KofamScan documented (F1 1.00 vs 0.93; FP 0 vs 14), edge
   concentrated at the shared-KO traps. ✅

---

## 7. Known limitations & residual errors

The training panel has zero residual FP/FN. The honest residuals live in the hold-out (12 FP +
15 FN), kept as evidence of where the model and curated biology diverge — **not** back-edited
into the GT:

1. **amoA_gamma ↔ methanotroph pmoA cross-fire** — *Methylococcus capsulatus* (the canonical
   amoA decoy) draws a γ-AOB HMM FP; the residual homology trap on novel γ-proteobacterial
   genomes. The single most informative hold-out finding.
2. **Nitrobacter nxrA/nxrB ≈ narG sequence-inseparability** — resolved by operon synteny on
   nucleotide/MAG input (§3.7); on bare proteomes (no gene coordinates) it persists. A divergent
   Chloroflexi NOB (*Nitrolancea*) also tripped the synteny resolver's narI-proximity heuristic
   (a narG FP) — the signature may not be reliable outside the Proteobacteria.
3. **Strict-GT FPs** — anammox NXR-like proteins (*Scalindua*), facultative nap+nar / nir
   redundancy (*Stutzerimonas stutzeri* F2a), methanotroph NO-detox NorB (*Methylococcus*):
   cases where the literature-conservative GT may itself be the stricter party.
4. **Annotation-gap FNs** — concentrated on *Anaeromyxobacter dehalogenans* (9 of 15); the
   unreviewed reference proteome appears to under-annotate variants the literature documents.

These are triageable to lift the hold-out F1 (per-cell literature review) but are documented
honestly rather than tuned away.

## 8. Conclusion
On a 31-genome training panel with curated, gene-content-verified, contamination-screened
ground truth, ncycle-pipeline achieves **micro-F1 1.00 [1.00, 1.00]** with **zero false
positives/negatives**, **perfect homology-trap precision (1.00)**, and all 9 pathways at
F1 = 1.00. On an 8-genome / 5-phyla independent hold-out it generalizes at **F1 = 0.85
[0.70, 0.95]** — a credible CI with real empirical width. It outperforms raw KofamScan
(F1 1.00 vs 0.93; FP 0 vs 14), with the advantage concentrated exactly at the shared-KO
homology traps the tool was designed to resolve. A CI-aware regression gate and full
provenance pinning protect these metrics going forward. The pipeline is in a
publication-ready state; remaining work-lines are elective (manuscript, external benchmark,
sister-pipeline parity) — see `ROADMAP.md` § Next-direction options.
</content>
</invoke>
