# ncycle-pipeline — Hardening Roadmap (MVP → robust)

> **Status (2026-06-10, Phase 5 COMPLETE; headline reframed after a statistics audit):**
>
> **Headline generalization: de-leaked HOLD-OUT F1 = 0.84 [0.67, 0.95]** on the
> 8-genome / 5-phyla hold-out (188 cells; α-/δ-/γ-proteo, Chloroflexi,
> Planctomycetes), after excluding 6 seed-contaminated positive cells
> (`detect_seed_leakage.py`). The seed-contaminated ("as-is") hold-out is
> 0.85 [0.70, 0.95] — leakage moves it only +0.01, so the contamination is
> real but small in effect. This is the credible real-world number.
>
> **Training F1 = 1.00 (0 errors / 645 cells)** is an **IN-SAMPLE fit** — the
> panel was used for threshold/seed calibration, so it reports training error,
> not generalization. The [1.00, 1.00] bootstrap CI is **degenerate** (a
> zero-error scope can't produce variance); the honest statement is "0 errors
> is consistent with a true per-cell error rate ≤ 0.47% (rule of three, 95%)."
> Do not headline the 1.00.
>
> **Stats-audit caveats (2026-06-10, see REPORT § Audit 2026-06-10):**
> - GT was tuned alongside detection on the training panel (~4/23 cell-flips
>   moved a label toward present in step with a new detector, all externally
>   documented); the rest were biology corrections, audit retractions, or new
>   true-negative cells. Real but bounded circularity — affects the 1.00 only.
> - Several per-target perfect scores rest on 1–2 positive cells (19/49 targets
>   flagged `f1_reliable=false`); "all traps solved" is carried by thin n_pos.
>
> **Acceptance criteria (reframed):**
> 1. Coverage — every pathway ≥1 training TP + every clade ≥1 hold-out TP.
> 2. Accuracy — F1 with 95% genome-cluster bootstrap CIs (B=10,000, seed=1234);
>    zero-error scopes additionally report a rule-of-three error bound.
>    Generalization carried by the de-leaked hold-out (0.84), not the 1.00.
> 3. Traps separated — amoA↔pmoA, nxr↔nar, nirB/D/A↔dsr clade-discriminated;
>    trap **precision** (robust to seed circularity — FPs land on never-seeded
>    GT-absent cells) is the trustworthy trap metric.
> 4. Benchmark — beats raw KofamScan + 3 tools on the training panel, 8/8
>    pre-registered BH-FDR contrasts; most credible component is trap precision.
>
> **Phase 5 deliverables** (in ~chronological order):
> - P5.0 — Burkholderia contamination of `test_panel/Nwinogradskyi_Nb255.faa`
>   resolved (file replaced with real Nb-255, GT bumped v10→v11 with v6 nosZ
>   + v9 norZ retractions, panel-QC gate `validate_panel.py` shipped).
> - P5.1 — bootstrap CIs, nosZ_clade2 TC hold-out leakage fix, GT-vs-detection
>   decomposition (`decompose_gt.py`), per-target F1 denominator capping,
>   UniProt seed snapshot cache + Pfam pin + `verify_seeds.py`.
> - P5.2 — Fragment AmoA seeds replaced, nrfH + hzsB BLAST seeds populated,
>   Chloroflexota nosZ_clade2 seeds added (TC 418 → 500), norB/norZ mutex
>   verified clean. **0/5 audit-supplied UniProt accessions were correct;
>   the seed cache caught every error at fetch time.**
> - P5.3 — synergy/complex any-of slots for clade paralogs, GT v12
>   programmatic TN expansion, hold-out 2→8 genomes spanning 5 phyla
>   (the headline CI-width fix).
> - P5.4 — BLAST gate identity-min documented, KO-filtered BLAST utility
>   (~50× cost cut), UniRef90 dead-code cleaned, CI-aware regression
>   floors (replaces point-estimate `MIN_HOLDOUT_F1`).
>
> **Pipeline is in publication-ready state.** Full narrative in
> `validation/REPORT.md` § Status + § Audit 2026-05-30 (sections
> "Post-audit response (P5.0)" through "Post-audit response (P5.4)").
> No active work-line; next-direction options listed at the foot of
> Phase 5 below.

Path from the working MVP (KO-primary detection, end-to-end on an 8-genome sanity
panel) to a defensible, measured tool. The MVP proved the architecture; this
roadmap closes the accuracy, coverage, and disambiguation gaps. Most machinery
already exists and is reused from `Holomicrobiome-ewaste/ewaste-pipeline`.

## Definition of "robust" (acceptance criteria)

1. **Coverage** — every one of the 9 pathways has ≥1 true-positive reference genome
   in the panel, and every target has ≥1 positive and ≥1 negative test genome.
2. **Accuracy** — per-target precision/recall/F1 with 95% bootstrap CIs; headline
   **median F1 ≥ 0.90**; no pathway's diagnostic marker below F1 0.80.
3. **Traps separated** — `nxrA`↔`narG`, `nxrB`↔`narH`, `amoA`↔`pmoA`,
   `nirB/D/A`↔`dsrAB` each resolve correctly on a panel containing *both* sides
   (precision ≥ 0.9 on the shared-KO/Pfam targets).
4. **Benchmark** — pipeline F1 ≥ raw KofamScan and ≥ NCBIfam on the same panel,
   with multiple-testing-corrected significance on the wins.

## Where we are (MVP baseline, measured by coverage)

- 32/49 targets confirmed in ≥1 genome; 6 pathways exercised with positives.
- **3 pathways unexercised** (no positive genome): nitrite oxidation (NOB),
  anammox, comammox; plus AOA (`amoA_archaeal`) and eukaryotic `NR`.
- **No accuracy metric exists** (no ground truth → no F1).
- Known error sources: BLAST seeds are auto-fetched by gene name (uncurated) →
  `amoC`/`napA`/`nirK` over-disqualified (false negatives); a divergent
  methanotroph's `amoA` single-gene false-confirms; `narrow-no-IPR` noise.

---

## Gap inventory → fix → reusable asset

| # | Gap (from the run) | Fix | Reuse |
|---|---|---|---|
| G1 | No accuracy metrics | Ground-truth panel + scoring harness | `validation/scripts/*` (all target-agnostic) |
| G2 | 3 pathways + AOA/NR never tested | Add NOB, anammox, comammox, AOA, a eukaryote, and decoys (methanotroph, SRB) to the panel | NCBI `datasets` (env present) |
| G3 | Auto-fetched BLAST seeds → false neg/pos | Curate + verify per-target UniProt seeds; tune `blast_identity_min` | `build_blast_db.py`; ewaste `targets/*/manifest.yaml` discipline |
| G4 | Shared-KO/Pfam traps not truly separated | Train clade-specific custom HMMs + calibrate TC | `build_custom_hmms.py`, `calibrate_tc.py` (already in repo) |
| G5 | KOfam `domain`-type thresholds compared to full score | Pass `score_type` through; use domain score for domain-type KOs | `apply_rules.py` `parse_hmmscan_domtbl` |
| G6 | `amoC` complex never reaches 100% even in true AOB | Falls out of G3 (seed fix) + G4 | — |
| G7 ✅ | Reuses `ewaste-pipeline` conda env; KOfam DB unpinned | **DONE (2026-05-25):** `run.py` defaults to the standalone `ncycle-pipeline` env (`NCYCLE_ENV` override); KOfam pinned by sha256 + verified at build | `envs/ncycle.yaml` |

---

## Phased plan (critical path: P1 → P2 → P3 loops)

### Panel expansion + hold-out validation ✅ DONE (2026-05-25)
Panel **16→33 genomes** (≥3/pathway + 3 negatives) via specialist audit; 2 independent
hold-outs (*Bradyrhizobium*, *Rhodopseudomonas* — genera absent from training).
**Training F1 0.95 (P 0.99); HOLD-OUT P 1.00, R 0.94, F1 0.97** → tool generalizes, not
overfit. Empirical verification caught 3 genome-identity bugs (japonica, M. extorquens
mislabel, Klebsiella nif-neg). See `validation/REPORT.md`.

### P1 — Validation harness + baseline ✅ DONE (2026-05-24)
**Baseline: ALL micro-F1 0.84; non-trap 0.92; homology-trap 0.54** (15-genome
panel, 273 curated cells). See `validation/REPORT.md`. The accuracy deficit is
concentrated exactly in the trap targets + auto-fetched seeds → confirms P2/P3
priorities. Harness reused: `validation/build_ground_truth.py`, `score_ncycle.py`.

Establishes the measuring stick before any tuning, exactly as ewaste did (B1 baseline).
1. **Panel** (G2): assemble ~40–60 reference proteomes covering all 9 pathways
   **plus** the missing positives — *Nitrobacter*/*Nitrospira moscoviensis* (NOB),
   *Nitrospira inopinata* (comammox), *Kuenenia*/*Brocadia* (anammox),
   *Nitrosopumilus* (AOA), a eukaryote (NR), and decoys (*Methylococcus*
   methanotroph for amoA/pmoA, *Desulfovibrio* SRB for nir/dsr). Keep the current 8.
2. **Ground truth** (G1): adapt `build_ground_truth.py` to N-cycle (UniProt/KEGG
   annotation → present/absent per genome) + a `LITERATURE_OVERRIDES` block for
   known gaps. Output `validation/ground_truth.tsv`.
3. **Score** (G1): run `score_phase1.py` + `bootstrap_ci.py` on the MVP calls →
   first per-target F1 table + CIs. Run `contamination_check.py` (panel/seed overlap).
4. **Acceptance:** a baseline metrics table; ranked list of weakest targets.

### P2 — Seed curation ✅ DONE (2026-05-25)
**ALL micro-F1 0.84 → 0.94; trap 0.54 → 0.80; FP 14 → 2.** Curated clade seeds for
nxrA/nxrB/napA/nirK/hdh/amoC/amoA + fixed the `ancient()` diamond-DB bug (DB rebuilds
now propagate to BLAST). Residual failures are now genuine P3 needs (comammox amo &
Nitrobacter nxr below-threshold / not BLAST-separable). See `validation/REPORT.md`.

Iterate per-target, re-scoring after each batch (ewaste's B-loop method):
1. Replace empty `blast_refs_uniprot` with curated, function-verified, clade-spread
   UniProt accessions (drop fragments/outliers); set `blast_identity_min` per target.
2. Prioritize the measured false-negatives (`amoC`, `napA`, `nirK`, `nrfH`) and the
   gated traps.
3. **Acceptance:** false-negative gates resolved; F1 improves vs P1 baseline on the
   affected targets; obligatory complexes reach 100% in true positives.

### P3 — Custom HMMs for the homology traps ✅ ROUND 1 DONE (2026-05-25)
**ALL F1 0.94→0.95; trap 0.80→0.86; ammonia-ox 0.71→0.80.** Built+calibrated nxrA/nxrB
(Nitrospira/Nitrobacter NOB) + amoA (AOB+comammox) HMMs; integrated (precedence
custom>KO>Pfam). Caught a panel data bug (the "comammox" genome was N. japonica NOB →
GT corrected to curated_v2). **Key finding:** Nitrospira-type NXR is cleanly separable
(score 2276 vs denitrifier ~530) but **Nitrobacter NxrA≈NarG (543 vs 544) is NOT
sequence-separable** — needs operon/phylogenetic context (documented limit; the sole
residual nxr FN + narG FP). **P3b: added a verified comammox *N. inopinata*** → amoA/nxrA/nxrB **confirmed via
custom-hmm** (KO alone misses them), comammox synergy complete. **All 9 pathways now
have a validated positive; ALL F1 0.95, trap 0.88, nitrite-ox 0.80, ammonia-ox 0.89.**
Remaining: AOA clade HMMs, nrfH/norB tweaks, Nitrobacter nxr (needs phylogeny). See `validation/REPORT.md`.

The shared-KO/shared-Pfam families that BLAST gating alone can't fully resolve:
1. Train clade HMMs with `build_custom_hmms.py` (CD-HIT cluster → MAFFT → hmmbuild)
   and calibrate cutoffs with `calibrate_tc.py` for: `nxrA`/`nxrB` vs `narG`/`narH`
   (NOB vs denitrifier clades), `amoA` AOB/AOA/comammox clades vs `pmoA`, and
   `nirB/nirD/nirA` vs `dsrAB`. Drop the trained `targets/{tid}/{tid}.hmm` in place —
   `build_hmm_db.py` already splices them and `apply_rules.py` already prefers
   custom > KO.
2. Fix G5 (domain-type KOfam thresholds) while in `apply_rules.py`.
3. **Acceptance:** criteria #3 (traps separated) met on a both-sides panel; the
   methanotroph `amoA` false-confirm eliminated.

### P4 — Benchmark + CV ✅ DONE (2026-05-25)
Fixes (nrfH ko_tc override→70; dropped harmful nirA/nirK gates) → **training F1 0.97
(P 0.99)**, DNRA 0.96, denitrification 0.97. **Hold-out P=R=F1=1.00.** **Comparator:
pipeline 0.97 vs raw KofamScan 0.93** (FP 3 vs 14; trap wins nxrA 0.75 vs 0.18, nxrB
0.75 vs 0.22, amoA 0.89 vs 0.73). **LOGO-CV:** nxrA generalizes across held-out
Nitrospira genus (985–996 ≫ TC700); amoA clade-sensitive (comammox seeds needed — not
overfit). Tools: `validation/compare_kofam.py`, `ko_tc` override in build_hmm_db.
Former residuals, both now addressed (2026-05-25): **Nitrobacter nxr↔narG resolved by
operon synteny** (`resolve_nxr_synteny` in apply_rules — narI/K00374 proximity; nucleotide/
MAG input only; validated on N. winogradskyi + E. coli). **AOA amoA_archaeal** fixed by
dropping the empty-seed BLAST gate (F1 0.00→1.00). See REPORT "Operon-synteny resolution".

### Phase 5 — Post-audit corrections + publication readiness (pick up here next session)

**Why this phase exists:** the 2026-05-30 dual-agent review (computational-biology
statistics + environmental-microbiology genomics specialists, audit reports
captured in `validation/REPORT.md` § Audit 2026-05-30) flagged a publication-blocking
contamination + several rigor / biology issues that need to clear before the
ALL F1 = 1.00 claim is defensible. Work the items in the order listed —
P5.0 is blocking, P5.1-P5.2 are publication-readiness, P5.3-P5.4 are
nice-to-have polish.

#### P5.0 — ✅ DONE (2026-05-30): Burkholderia contamination + Nwinogradskyi re-verification

**Resolution:** file replaced with real Nb-255 (NC_007406.1, 3262 proteins) →
pipeline re-run on full panel → GT bumped v10 → v11 (v4 narG STANDS, v6 nosZ
RETRACTED, v9 norZ RETRACTED, nosZ_clade2 → A added) → re-scored. **ALL F1
= 1.00 verified** against v11 GT on the corrected panel; 0 residual FN/FP.
See `validation/REPORT.md` § Audit 2026-05-30 → "Post-audit response (v11)".
Panel-QC gate `validation/validate_panel.py` landed (Makefile target
`validate-panel`, wired as a `regression` / `regression-score` dependency) —
fails fast on a Burkholderia-style swap (verified against the quarantined
contaminated file: 42/42 wrong-genus → exit 1). Historical context below.

**The fact:** `test_panel/Nwinogradskyi_Nb255.faa` is *Burkholderia thailandensis*
(5607 proteins, 5258 explicitly `[Burkholderia thailandensis]`-tagged, ZERO
Nitrobacter tags). Verified independently.

1. **Replace the file.** Pick one of:
   - Re-prodigal `test_synteny/Nwinogradskyi.fna` (NC_007406.1, the real Nb-255
     genome; existing prodigal output is at `results/Nwinogradskyi/prodigal/Nwinogradskyi.faa`
     = 3262 proteins). This is the immediate-fix path — proteome already exists.
   - Download RefSeq GCF_000012685.1 protein FAA for the canonical RefSeq
     annotation (~3120 proteins). Cleaner provenance for publication.
2. **Re-run the pipeline on the corrected Nb-255 row only**
   (`python run.py --input ../test_panel --skip-db-setup`).
3. **Re-evaluate the three Nb-255 GT corrections** against the real proteome:
   - **v4 narG A→P** — literature (Starkenburg 2006/2008) supports a narGHI operon
     in real Nb-255 at 0.86 Mb on NC_007406.1; likely to hold but verify via the
     operon-synteny output.
   - **v6 nosZ A→P** — Starkenburg supports facultative denitrification but the
     specific evidence used (WP_080511513.1 = Burkholderia TAT-dep nosZ) was wrong;
     re-verify with real Nb-255 K00376 + PF18764 hits.
   - **v9 norZ A→P (qNor)** — the "new finding" evidence (WP_009890104.1) is
     Burkholderia thailandensis qNor. Real Nb-255 is not a documented qNor carrier;
     **expect to retract this correction**. Restore the v8 docstring framing.
4. **Bump GT to v11** with explicit retract-or-confirm comments on each Nb-255
   line. Update `Nwinogradskyi_Nb255` entry; update REPORT.md § Batch 3c
   "Nwinogradskyi qNor finding" with a retraction footnote if v9 norZ flips back.
5. **Re-run + re-score the panel.** The ALL F1 value will move; document the
   new value with honest provenance (separate "GT corrections" vs "detection
   improvements" contributions — see P5.1.3).
6. **Add panel-QC gate** (`validation/validate_panel.py`): per-genome, sample 50
   protein headers and fail if <60% match the expected genus/species from the
   filename. Wire into the Makefile + run before any GT regeneration.

Acceptance: real Nb-255 proteome installed; v4/v6/v9 cells re-evaluated with
documented outcomes; panel-QC gate in place to prevent recurrence; updated
REPORT.md + ROADMAP.md banner with the new (not 1.00) headline F1.

#### P5.1 — Publication-readiness statistical rigor

Audit-flagged: acceptance criterion #2 ("95% bootstrap CIs on F1") is **literally
unmet** in the headline path. Machinery exists but isn't wired in.

1. ✅ **DONE (2026-05-30).** Wired genome-cluster bootstrap CIs into
   `validation/score_ncycle.py` (B=10,000, seed=1234, percentile method
   mirroring `validation/benchmark/benchmark_stats.py:130-144`). CIs print
   inline with every aggregate, per-pathway, and hold-out point estimate.
   All CIs collapse to [1.00, 1.00] on curated_v11 because the panel has
   zero FP + zero FN globally (honest structural behavior); synthetic-
   imperfect-panel sanity check confirms the function produces real width
   on imperfect data. Hold-out CI remains structurally uninformative
   (n=2 genomes); the meaningful width question is deferred to P5.3.5
   (≥8-genome hold-out).
2. ✅ **DONE (2026-05-30).** Recalibrated `nosZ_clade2` TC 420 → 418
   excluding hold-out genomes. New ceiling: Smeliloti 414.2 (training-only) +
   3.8 margin. Hold-outs (Bdiazoefficiens 416.8, Rpalustris 414.6) excluded
   from calibration logic but both still correctly classified absent at
   TC=418. No panel classification change; ALL F1 = 1.00 [1.00, 1.00],
   hold-out F1 = 1.00 [1.00, 1.00] preserved. Manifest updated
   (`targets/nosZ_clade2/manifest.yaml`), HMM DB rebuilt, panel re-run +
   re-scored. v10 narrative cited a Nwinogradskyi 414.5 figure that was a
   contamination artifact (Burkholderia cross-fire); real Nb-255 has zero
   `nosZ_clade2` hits at any E-value. See REPORT § Audit 2026-05-30 →
   "Post-audit response (P5.1.2)".
3. ✅ **DONE (2026-05-30).** Built `validation/decompose_gt.py` which
   reconstructs each historical GT (v3 → v11) by walking documented
   cell-flips in reverse and re-scores the current call matrix against
   each. **Hold-detection-constant trajectory:** v3 GT F1 = 0.984, v11 GT
   F1 = 1.000 → GT revisions account for +0.016 F1 points. Most of the
   improvement from P1 baseline (F1 ≈ 0.84) is therefore detection-side.
   **Hold-GT-constant leg is not retroactively computable** without baseline
   pipeline-state snapshots — the chronological F1 trajectory in REPORT.md
   sections P1 → Batch 3d serves as observational record for that leg.
   Attribution: 5 biology corrections + 13 detection-enabled (3 prediction
   wins + 10 TN-coverage expansions from new targets, the latter don't
   affect F1) + 5 audit retractions. See REPORT § Audit 2026-05-30 →
   "Post-audit response (P5.1.3)".
4. ✅ **DONE (2026-05-30).** `MIN_RELIABLE_POSITIVES = 3` floor wired into
   `validation/score_ncycle.py`. Per-target rows now carry `n_pos` (TP+FN)
   and `f1_reliable` (boolean) fields. Stdout flags low-reliability rows
   with `*` and prints a legend listing the 8 affected targets; TSV gets
   `n_pos` and `f1_reliable` columns. Aggregates unchanged. Regression
   gate still PASS. The two new clade HMMs (amoA_gamma, nosZ_clade2) are
   in the flagged list — biological positives would need to be added to
   the panel (P5.3.5 hold-out expansion) to clear them. See REPORT § Audit
   2026-05-30 → "Post-audit response (P5.1.4)".
5. ✅ **DONE (2026-05-30).** UniProt seed snapshot cache + Pfam-A sha256
   pin shipped. New: `workflow/scripts/_seed_cache.py`,
   `workflow/scripts/bootstrap_seed_cache.py` (one-shot pin from current
   build), `validation/verify_seeds.py` (read-only drift check;
   `make verify-seeds`). Modified: `build_blast_db.py` (cache-first with
   `--refresh-seeds`), `build_hmm_db.py` (Pfam pin verification on fetch).
   Pinned: 124 unique UniProt accessions across 32 targets + PF12942
   (Pfam-A archaeal amoA fallback). UniRef90 path is dead code per
   P5.4.3 — out of scope. Drift contract: silent UniProt revisions are
   now caught by `make verify-seeds` before they can invalidate TC
   calibration. See REPORT § Audit 2026-05-30 → "Post-audit response
   (P5.1.5)".

Acceptance: headline F1 reported with CIs; nosZ_clade2 TC clean of hold-out
leakage; GT-vs-detection table present in REPORT methods; per-target
denominators visible; full provenance pinning.

#### P5.2 — Biology seed-quality fixes (env-microbio audit findings)

Audit-flagged issues with the canonical seed sets — each is a one-line
config edit + an HMM/BLAST DB rebuild.

1. ✅ **DONE (2026-05-30).** Replaced Q7WWN6 + A0A1P8VZN2 (both N. multiformis
   Fragments) with Q51142 (N. multiformis AmoA1, 274 aa full-length) +
   A0A0F7KK52 (N. communis, 274 aa full-length). amoA HMM re-trained
   (5 candidates → CD-HIT → 4 → MAFFT → hmmbuild); TC=400 stands
   (positives 511–590, cross-fire 207–280, gap unchanged). Two audit-error
   findings surfaced during verification: **Q7WWN7** (audit-suggested
   N. eutropha replacement) was *another* fragment; **D3RUW6** (audit-
   suggested N. briensis AmoA) was ATP synthase from Allochromatium vinosum
   (wrong protein + organism). Both dropped from the swap; documented in
   REPORT + manifest. Final 5-seed BLAST set all full-length. F1 = 1.00
   [1.00, 1.00] unchanged. See REPORT § Audit 2026-05-30 → "Post-audit
   response (P5.2.1)".
2. ✅ **DONE (2026-05-30).** `nrfH.blast_refs_uniprot` populated with 4
   verified reviewed seeds: Q9S1E6 (Wsuccinogenes NrfH), P0ABL1 (E. coli
   K12 NrfB), Q72EF4 (D. vulgaris Hildenborough NrfH), P45016 (H. influenzae
   NrfB). **Both audit-named accessions were wrong:** P0C278 = Shewanella
   fumarate reductase (not Wsuccinogenes NrfH); P0AAJ7 = Shigella formate
   dehydrogenase (not E. coli NrfH); Shewanella oneidensis NrfH doesn't
   exist in UniProt at all. Audit's organism guidance was correct; only
   the accessions were wrong. **Pre-existing state correction:** nrfH was
   already F1=1.00 via KO K15876 at ko_tc=70 (NOT silently dead — the
   audit's `requires_blast_for_confirmation: true` claim was outdated; the
   target has `blast_fallback: true` only). The seed population adds an
   active confirmation layer that was previously inert. BLAST DB rebuilt
   (80 → 84 seeds); pipeline + regression unchanged at F1=1.00 [1.00, 1.00].
   See REPORT § Audit 2026-05-30 → "Post-audit response (P5.2.2)".
3. ✅ **DONE (2026-05-30).** `hzsB.blast_refs_uniprot` populated with 3
   verified anammox seeds: Q1Q0T4 (Kuenenia stuttgartiensis HzsB, reviewed
   sp/), A0A0M2V157 (Ca. Brocadia fulgida HzsB), A0A1E3XDZ9 (Ca. Scalindua
   rubra HzsB). Jettenia only exists as 128-aa PCR fragment; other
   "Hydrazine synthase β" UniProt hits from Xanthomonas/Pseudomonas/
   Burkholderia are annotation errors on non-anammox bacteria. F1=1.00
   unchanged. See REPORT § Audit 2026-05-30 → "Post-audit response (P5.2.3–
   P5.2.5)".
4. ✅ **DONE (2026-05-30).** Added 2 Chloroflexota Sec-dependent nosZ seeds
   to `targets/nosZ_clade2/`: A0A7C1FJ72 (Caldilinea aerophila, Caldilineae
   class type strain) + A0A7C1K3A7 (Thermomicrobium roseum, Thermomicrobia
   class type strain). HMM re-trained (5 → 7 seeds). Cross-fire scores
   shifted up 30-50 bits uniformly → **TC recalibrated 418 → 500**
   (training-only ceiling 456.9 + 2.3σ; no hold-out leakage; ~590-bit
   discrimination gap to clade-II positives preserved). 6 phyla covered.
   F1=1.00 unchanged. See REPORT § Audit 2026-05-30 → "Post-audit response
   (P5.2.3–P5.2.5)".
5. ✅ **DONE (2026-05-30).** Inspected `results/multisample_matrix.tsv`
   (33 genomes): **zero cells have both norB and norZ at predicted-present
   status.** The existing BLAST gates already provide mutual exclusion at
   the sequence-identity level (norB's β-proteo cNor seeds at identity≥55,
   norZ's qNor seeds at identity≥40 — sequence-distinct enough that a
   single K04561 protein hits one and is BLAST-disqualified by the other).
   Empirically clean: Cnecator + Synechocystis → norZ=confirmed +
   norB=disqualified; other denitrifiers → norB=confirmed + norZ=disqualified.
   **No explicit mutex code needed.** See REPORT § Audit 2026-05-30 →
   "Post-audit response (P5.2.3–P5.2.5)".

Acceptance: no Fragment seeds anywhere; all "requires_blast" targets have
populated `blast_refs_uniprot`; nosZ_clade2 has Chloroflexi coverage; norB/norZ
mutex verified.

#### P5.3 — Capability extensions (long-flagged, mostly Batch 3 follow-ups)

Items previously deferred from Batch 3 follow-ups (REPORT § Block 3 follow-ups).

1. ✅ **DONE (2026-05-30).** `nitrifier_denitrification` synergy updated to
   `requires: [[amoA, amoA_gamma, amoA_archaeal], [nirK, nirS], [norB, norZ]]`.
   Noceani's γ-AOB nitrifier-denit phenotype now reports 1.000 complete via
   amoA_gamma+nirK+norB. See REPORT § Audit 2026-05-30 → "Post-audit response
   (P5.3.1–P5.3.3)".
2. ✅ **DONE (2026-05-30).** 4 nosZ-bearing synergies updated + 1 refined:
   complete_denitrification + n2o_sink + n2o_emitter (forbids) + n2o_sink_only
   accept `[nosZ, nosZ_clade2]` any-of; nosZ_clade_I_likely kept clade-I-
   specific but gained `forbids: [nosZ_clade2]` (phylogeny-rigorous now that
   the clade-II HMM exists). nitrate_to_nitrite_leak doesn't involve nosZ —
   no change. See REPORT.
3. ✅ **DONE (2026-05-30).** Chose option (b): extended
   `compute_complex_completeness.py:75-90` to support any-of slots for complex
   members (mirrors the synergy evaluator pattern). Updated
   `ammonia_monooxygenase` complex `members: [[amoA, amoA_gamma, amoA_archaeal],
   amoB, amoC]`. Also flatten-on-read fixes in `make_pathway_heatmap.py` +
   `make_applications_panel.py` so downstream consumers handle list members.
   Noceani complex now 1.000 complete via amoA_gamma+amoB+amoC. F1=1.00
   unchanged; regression 15/15 PASS. See REPORT.
4. ✅ **DONE (2026-05-30).** Programmatic TN-expansion post-pass added to
   `validation/build_ground_truth.py` (`EXPAND_TN_TARGETS` + source tag
   `expanded_tn_v12`). amoA_gamma cells 6 → 33; nosZ_clade2 cells 5 → 33.
   Total GT 645 → 700 cells. F1=1.00 unchanged. P5.1.4 reliability flag
   not yet cleared (depends on n_pos ≥ 3; biological positives come from
   P5.3.5 hold-out expansion). See REPORT § Audit 2026-05-30 → "Post-audit
   response (P5.3.4)".
5. ✅ **DONE (2026-05-30).** Hold-out expanded 2 → 8 genomes spanning 5 phyla
   (audit minimum ≥4). New: Adehalogenans 2CP-1 (δ-proteo, clade-II NosZ),
   Nhalophilus Nc4 (γ-AOB), Nhollandica Lb (Chloroflexi NOB), Sbrodae
   (Planctomycetes anammox), Sstutzeri F2a (γ-proteo full denit),
   Mcapsulatus Bath (γ-proteo methanotroph amoA decoy). **Hold-out F1 =
   0.85 [0.70, 0.95]** — first real CI width. amoA_gamma + nosZ_clade2
   now have hold-out positives. Training F1 = 1.00 [1.00, 1.00] preserved.
   12 FP + 15 FN documented as honest detection/annotation findings (NOT
   back-edited into GT, which would defeat hold-out purpose). Regression
   floor `MIN_HOLDOUT_F1` lowered 0.95 → 0.70 to match new CI lower bound;
   P5.4.4 will switch to CI-aware floor. See REPORT § Audit 2026-05-30 →
   "Post-audit response (P5.3.5)".

Acceptance: synergies and complexes accept all clade variants; GT TN coverage
panel-wide for the two new clade HMMs; hold-out has ≥8 genomes with
documented CIs.

#### P5.4 — Efficiency + reproducibility hardening (lowest priority)

Audit-flagged items that aren't blocking but improve scalability /
reproducibility for the publication and downstream use.

1. ✅ **DONE (2026-05-31).** Banner + 8 per-target `gate_rationale:` fields
   added. Reviewers can read why every gate value sits where it does. See
   REPORT § Audit 2026-05-30 → "Post-audit response (P5.4)".
2. ✅ **DONE (2026-05-31); WIRED + RE-EVALUATED (2026-06-09).**
   `workflow/scripts/filter_proteome_by_hmm.py` shipped as optional utility
   (retention ~1.57% over the 39-genome panel → ~63× fewer DIAMOND queries).
   Now wired into the DAG behind config `options.filter_blast_query_by_hmm`
   (rule `filter_proteome_by_hmm` → both `diamond_blastp_*` rules query the
   filtered FASTA; hmmscan + apply_rules synteny keep the FULL proteome).
   **Toggle DEFAULTS OFF — the filter is NOT prediction-neutral here:** the
   BLAST gate is complementary to the HMMs, so BLAST-only targets (assimilatory
   `narB`/K00367, `nasD`/K17877) detected on divergent proteins with zero
   hmmscan hits get dropped → training F1 1.00 → 0.996 (2 FN, narB E. coli +
   nasD N. europaea), regression FAILS when ON. With the default OFF the panel
   re-runs bit-identical (matrix unchanged, 15/15 regression PASS). Lesson: the
   filter only speeds the already-cheap DIAMOND-vs-seeds step, not the hmmscan
   bottleneck that dominates large runs — so it is **not** the GTDB-scale lever;
   left wired + toggleable for throughput experiments only. The real scaling
   lever is hmmscan throughput (per-genome parallelism / candidate prefilter).
3. ✅ **DONE (2026-05-31).** All UniRef90 dead code removed from
   `build_custom_hmms.py` (`cmd_expand`, `ensure_uniref90`,
   `_extract_records_by_acc`, `_organism_from_header`,
   `_parse_phmmer_tblout`, UniRef90 constants). Script size 450 → 228 lines.
   `cmd_build` auto-fetches missing accessions from UniProt at build time,
   making `expand` structurally unnecessary. Smoke-tested; regression PASS.
4. ✅ **DONE (2026-05-31).** `validation/test_regression.py` switched to
   CI-aware floors (CI lower bound vs pinned v12 baseline + explicit slack).
   ALL F1 CI lo ≥ 1.00, hold-out F1 CI lo ≥ 0.65 (v12 observed 0.695 +
   0.045 slack), trap precision CI lo ≥ 0.95, per-pathway F1 CI lo ≥ 0.70.
   Output now shows `metric [lo, hi]` inline. Regression 15/15 PASS at v12
   baseline. P5.4.4 closes acceptance criterion #2.

Acceptance: BLAST gate thresholds documented; dead code paths cleaned;
regression test uses CI floors.

---

### Next-direction options (Phase 5 complete; no default work-line)

Pick one based on grant timing and intent. Each has a different scope
and audience:

1. **Manuscript draft.** The methods narrative is now defensible end-to-
   end. Sections to write: introduction (why a curated N-cycle pipeline
   over raw KofamScan), methods (panel composition + curation + bootstrap
   CI methodology + clade-discrimination layers), results (training F1,
   hold-out CI, comparator panel), discussion (limitations: amoA decoy
   class on novel genomes per P5.3.5 FPs; future directions). Estimated
   ~2-3 weeks for a draft.

2. **External benchmark at scale.** P5.4.2 made the ~50× DIAMOND cost
   cut possible. Sweep against MetabolicHMM, NCycDB, DRAM-v on 500-1000
   genomes (e.g. a GTDB representative slice). Would substantiate the
   "vs raw KofamScan" claim with multiple comparators on a larger panel.
   Requires writing comparator adapters + a panel-curation script for
   GTDB selection. Estimated ~3-4 weeks.

3. **Ewaste-pipeline parity sweep.** Apply Phase 5 hardening to the
   sister ewaste-pipeline. Most of P5 is target-agnostic
   (seed cache + verify_seeds + decompose_gt + CI-aware regression +
   filter_proteome_by_hmm + validate_panel). The non-portable parts
   are the clade-specific HMMs and per-target gate_rationale. Estimated
   ~1-2 weeks if ewaste's panel + GT structure mirrors ncycle's.

4. **Fresh sister-pipeline (carbon-cycle).** Build a new
   ccycle-pipeline from scratch using the now-established
   ncycle-pipeline template. Phase 5 hardening would be baked in from
   day 1. Sulfur is already done (scycle-pipeline). Estimated ~6-8 weeks
   for a comparable maturity level.

If none of the above resonates, possible follow-ups within ncycle:
- Investigate the 12 hold-out FPs + 15 FNs documented in P5.3.5 with
  per-cell literature triage (some may be genuine GT corrections, others
  detection limits worth fixing). Could lift hold-out F1.
- Tighten amoA_gamma vs methanotroph pmoA discrimination — the Mcapsulatus
  cross-fire P5.3.5 surfaced is a real residual trap.
- Add a Snakemake wire-in for `filter_proteome_by_hmm.py` to make
  P5.4.2 the default at panel sizes >100 genomes.

---

### P4 — original polish notes (G7)
1. `score_comparator.py` (vs NCBIfam) + `score_kofam.py` (vs raw KofamScan) +
   `correct_pvalues.py`; `logo_cv.py` leave-one-genus-out for the custom HMMs.
2. ✅ Standalone `ncycle-pipeline` conda env (`run.py` default, `NCYCLE_ENV` to override);
   KOfam release pinned by sha256 (`build_hmm_db.py` verifies the tarball + writes
   `resources/.cache/kofam_release.txt`; baseline = 2026-05-24, `b03d20b9…`).
3. ✅ Refine process-completeness logic (comammox, nosZ clade I/II, partial denitrification,
   nitrifier-denitrification) against the expanded panel. **DONE (2026-05-28):** synergy
   layer extended with `requires:`-as-slots (paralog alternatives — napA↔narG, nirK↔nirS,
   norB↔norZ, …), `forbids:` for named partials (n2o_emitter, n2o_sink_only,
   nitrate_to_nitrite_leak), `single_organism: true` to collapse misleading partials on
   non-community-completable phenotypes (comammox, anammox), and a heuristic
   `nosZ_clade_I_likely`/`n2o_sink_only` clade pair. Per-target accuracy unchanged
   (ALL F1 0.97, hold-out 1.00). Phylogeny-rigorous clade-II nosZ HMM deferred (future P5).
   See REPORT "Process-completeness refinement".
4. **Acceptance:** all four "robust" criteria met; a `validation/REPORT.md`.

---

## Sequencing notes & risks

- **Order matters:** P1 first (measure), then P2/P3 as scored iteration loops — never
  tune blind. P2 before P3 (cheaper wins first; some "traps" may resolve with seeds alone).
- **Mostly data, not code:** the harness, custom-HMM builder, and calibrator already
  exist and are target-agnostic — effort is curating panels/seeds and training ~10 HMMs.
- **Network:** panel + KOfam downloads are the throttled bottleneck (use the 12-way
  segmented download used for KOfam profiles).
- **Risk — anammox/Tier-3 multiheme** (`hzs`, `hdh`, `hao`): few reference genomes and
  no clean Pfam; may stay BLAST-only with modest F1 until more references exist
  (the lanmodulin precedent in ewaste). **Resolved (2026-05-25):** built+calibrated custom
  HMMs for the trapped octaheme markers `hao` (TC 622.8) and `hdh` (TC 741.6) — both fire as
  custom-hmm evidence. The `hzsA/B/C` subunits stay Tier-3 BLAST-only (per-subunit HMMs
  overlap each other + a KEGG↔UniProt subunit-naming conflict); they already score F1 1.00.
  `narG/narH` custom_hmm flags dropped (NxrA-inseparable). All F1 unchanged at 0.97.
  See `validation/REPORT.md` → "Tier-3 multiheme custom HMMs".
- **Risk — `amoB/amoC` genuinely shared** with methanotroph `pmoB/pmoC` at the KO level;
  full separation may require the complex-level call (amoA-clade gates the complex)
  rather than per-subunit certainty.
