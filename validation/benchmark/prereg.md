# Benchmark pre-registration — ncycle-pipeline vs external N-cycle tools

Fixed **before** running the external comparators, to keep the comparison honest
(no tuning either side after seeing results, no metric/endpoint shopping).

## Hypotheses
- **Primary (confirmatory):** on the homology-trap targets, ncycle-pipeline has higher
  **precision** than each comparator. Expected large effect (trap precision 1.00 vs the
  shared-KO failure mode ~0.2–0.5) → the panel is well powered here.
- **Secondary:** ncycle-pipeline has higher **overall micro-F1** than each comparator.
  Small expected effect (≈0.97 vs 0.93) → under-powered at n=33; reported, not relied on.
- **Exploratory:** per-pathway / per-step F1 (report CIs; not part of the corrected family).

## Comparators (default settings, no asymmetric tuning)
| tool | approach | resolves subunits? |
|---|---|---|
| **ncycle-pipeline** (this tool) | KO + custom HMM + BLAST gating (+ synteny on nt input) | yes (49 targets) |
| **raw KofamScan** | KO-only, stock `ko_list` thresholds, shared KOs → all targets | yes (data already available) |
| **METABOLIC** | HMM-based N/S/C cycle step detection | partial |
| **NCycDB** | curated N-cycle gene-family DB + DIAMOND | yes (gene families incl. nxrA/narG) |
| **DRAM** | MAG metabolism distillation (KEGG/Pfam) | coarse (steps/modules) |

Each tool's exact version + command line is recorded in `adapters.py` (and must be
filled in when run). Default thresholds only.

## Inputs
Same 33-genome panel proteomes for all tools. The 2 hold-out genera are **excluded**
from the cross-tool comparison (training panel only, as in `compare_kofam.py`).
Note: ncycle's operon-synteny step is dormant on proteome input — this benchmarks the
KO+HMM+gating core, which is the fair common ground.

## Comparison vocabulary (two pre-specified resolutions)
1. **Subunit** — the 49 ncycle targets. Coarse tools are mapped *down*: a step-level
   call is expanded to all member subunits (credit and blame applied symmetrically).
   The trap claim lives here.
2. **Step** — ~13 N-cycle transformation steps (`STEP_MAP` in `adapters.py`). ncycle is
   mapped *up* (step present iff its diagnostic marker present). Fairest to coarse tools.
   Reporting both guards against penalizing a tool merely for granularity.

Ground truth (`curated_v4`) is collapsed to each resolution by the same rules for all tools.

## Scoring & statistics
- **Scorer:** shared `compute_metrics()` logic; a cell counts only if it is in the
  ground truth. Positive class = "present".
- **Unit of independence = genome** (cells within a genome are correlated). All CIs use a
  **genome-level cluster bootstrap** (resample the 33 genomes with replacement, B=10,000).
- **95% CIs** (percentile) on precision / recall / F1 per tool, per subset, per resolution.
- **Primary significance:** paired **bootstrap of the difference** Δ = metric(ncycle) −
  metric(comparator), recomputed on the same resampled genomes each iteration; the win is
  significant if the 95% CI of Δ excludes 0 (two-sided bootstrap p reported).
- **Secondary significance:** **exact McNemar** on the cell-level discordant pairs
  (ncycle-correct/comparator-wrong vs vice-versa). Flagged secondary because it ignores
  genome clustering.
- **Multiple testing:** **Benjamini–Hochberg FDR** across the pre-specified family =
  {comparators} × {trap-precision (subunit), overall-F1}. Exploratory metrics not corrected.
- **Effect sizes** (ΔF1, Δprecision with CIs) reported alongside every p-value.

## Decision rules
- A comparator is "beaten on the trap" iff Δ trap-precision CI excludes 0 (BH q<0.05) AND
  the point estimate favors ncycle.
- Overall-F1 differences are reported with CIs; if a CI includes 0 we state the panel is
  underpowered for that contrast rather than claiming parity or a win.

## Known limitations
- n=33 genomes: well powered for the trap (large effect), modestly powered for the small
  overall gap. Panel expansion is the lever to power the secondary endpoint.
- Vocabulary mapping is a source of bias; it is fixed here, applied identically to all
  tools, and stored as a reviewable table.
- KEGG-based tools may have seen some panel genomes in their training references — a shared
  limitation that does not affect the trap-discrimination claim.

## Amendments (2026-06-07) — comparator set + reporting harmonization with scycle-pipeline
Recorded as a transparent post-hoc amendment; the **confirmatory family is unchanged** (trap
precision + overall-F1, BH-corrected, as above). These additions are descriptive only.
- **Comparator set finalized = {raw KofamScan, METABOLIC v4.0, DRAM v1.4.6}.** METABOLIC + DRAM
  were run (provenance in `adapters.py`); **NCycDB is DEFERRED to future work** (not installed),
  so the ncycle and scycle benchmarks share an identical comparator set. The `load_ncycdb`
  scaffold + `NCYC_MAP` are retained for a future run.
- **Three reporting frames, consistent with scycle-pipeline:** (A) whole-panel micro-F1;
  (B) per-(genome,target) independent-subset trap precision (`validation/trap_independence.py`);
  (C) train/hold-out split + leave-one-clade-out CV (`validation/logo_cv.py`). Frames A/B/C are
  descriptive cuts, NOT new confirmatory endpoints — the BH family is untouched. The cross-tool
  benchmark genome set is unchanged (training panel; hold-outs excluded as pre-registered).

## Amendment (2026-06-09) — NCycDB run (reverses the 2026-06-07 deferral)
The 2026-06-07 deferral of NCycDB is **withdrawn**: NCycDB (Tu et al. 2019) was installed and
run on the training panel. Provenance + the fixed family→target mapping are in `adapters.py`
(loader `load_ncycdb`, `NCYC_MAP`); the search is faithful to `NCycProfiler.PL` protein default
(`diamond blastp -k 1 -e 0.0001` vs `NCyc_100`). **The confirmatory family is unchanged** — NCycDB
simply adds its two pre-specified rows ({trap-precision, ALL-F1}) to the same BH-corrected family,
exactly as scaffolded; no endpoints, decision rules, genome set, or vocabulary rules were altered.
The `amoA_A`/`amoA_B` mapping was fixed from NCyc's own clade semantics, verified empirically
(amoA_A hit only by AOA archaea, amoA_B only by bacterial AOB/comammox/γ-AOB). Comparator set for
ncycle is therefore {raw KofamScan, NCycDB}; METABOLIC + DRAM remain available to slot in later via
the same harness. (This intentionally diverges from scycle's {kofam, METABOLIC, DRAM} set — noted
for transparency; the per-tool contrasts are independent, so the set difference does not affect any
single comparator's result.)

## Amendment (2026-06-09, later) — METABOLIC + DRAM run; comparator set complete
METABOLIC v4.0 and DRAM v1.4.6 were run on the training panel (provenance in `adapters.py`),
completing the comparator set to {raw KofamScan, METABOLIC, NCycDB, DRAM}. Both add their two
pre-specified rows ({trap-precision, ALL-F1}) to the same BH-corrected family; no endpoints,
decision rules, genome set, or vocabulary rules changed. All 8 contrasts (4 tools × 2 endpoints)
favor ncycle after BH-FDR (q ≤ 0.0012). METABOLIC + DRAM map through KEGG KOs / KEGG modules, so
their mapping tables are the existing KO→target map and `DRAM_STEP_MAP` respectively (fixed, not
panel-tuned). This supersedes the earlier 2026-06-07 deferral note and matches/extends scycle's
comparator set.

## Amendment (2026-06-14, B8) — battery harmonization + raw-KofamScan parser fix
Mirror of the scycle B8 amendment; **confirmatory family, decision rules, genome set, GT, and
vocabulary are unchanged** — reproducibility/methods fixes only.
- **Comparator-set parity:** ncycle = {KofamScan, METABOLIC, DRAM, NCycDB} and scycle =
  {KofamScan, METABOLIC, DRAM, SCycDB} — identical four-comparator shape with a domain DB on each.
  The 2026-06-09 note about ncycle "intentionally diverging" from scycle's comparator set is now
  **resolved**: scycle added SCycDB (its domain DB), so the two sister benchmarks are structurally
  identical and report the same battery the same way.
- **raw-KofamScan baseline parser fix:** the pipeline now emits `hmmsearch --domtblout` (KO in
  column 4, full-seq score in column 8); `load_kofam` had keyed on the legacy `hmmscan --tblout`
  column 1 and was silently scoring KofamScan all-absent. The adapter now detects either layout.
  Fix verified: the re-run reproduces the previously-reported KofamScan numbers exactly (ALL-F1
  0.916 [0.887,0.941], trap-precision 0.836) — confirming no drift in the other tools or in ncycle.
- **Final re-run:** B=10,000 genome cluster bootstrap, seed 1234, GT=`ground_truth.tsv`, n=31 panel
  genomes (hold-outs excluded) → `benchmark_results.tsv` + `COMPARISON_REPORT.md`. All 8 contrasts
  (4 comparators × {ALL-F1, trap-precision}) favour ncycle after BH-FDR (q<0.0012).

## Amendment (2026-06-16) — secondary TN-inclusive metrics (specificity + MCC)
Mirror of the scycle 2026-06-16 amendment. Post-hoc **secondary** analysis; the pre-registered
confirmatory family (trap-precision + ALL-F1, BH-FDR) is **unchanged** and bit-identical on re-run.
Added **specificity** = TN/(TN+FP) and **MCC** (all four confusion cells; robust to class imbalance),
which credit true negatives that F1 ignores. Same B=10,000 cluster bootstrap; the MCC family (ALL +
trap × 4 comparators) is BH-FDR-corrected **separately** from the primary family. All 8 favour ncycle
(q<0.0001). ncycle = 1.000 on both specificity and MCC.
