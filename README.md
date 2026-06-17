# ncycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **nitrogen cycle**.
Modeled on `Holomicrobiome-ewaste/ewaste-pipeline`; detection is **KO-primary**
(KOfam HMMs + adaptive per-KO thresholds) with custom clade HMMs and DIAMOND-BLAST
gating for the homology traps. Covers **49 targets** (47 KEGG Orthologies), **9
obligatory complexes**, and **7 process-completeness synergies** across all 9
N-cycle pathways.

**Status: validated.** End-to-end on the reference panel (training genomes + 8 independent
hold-out genera). The honest generalization headline is the **de-leaked hold-out micro-F1 ≈
0.84** (95% genome-cluster bootstrap CI [0.68, 0.95]; P 0.86, R 0.82). The **training-panel
F1 = 1.00 is an in-sample calibration fit** (thresholds/seeds tuned on those genomes; honest
per-cell error bound ≤0.47%, rule of three), **not** a generalization number. The
leakage-robust comparator advantage over raw KofamScan is **homology-trap precision**. See
[`validation/REPORT.md`](validation/REPORT.md) (§ Audit 2026-06-10 reframed the headline) and
[`ROADMAP.md`](ROADMAP.md).

## Run

```bash
cd ncycle-pipeline
python run.py --input <dir-of-.faa-or-.fna> --cores 8
# First run creates the standalone `ncycle-pipeline` conda env from envs/ncycle.yaml.
# To reuse an existing compatible env instead (e.g. the legacy ewaste one):
#   NCYCLE_ENV=ewaste-pipeline python run.py --input <dir> --skip-db-setup
```

`run.py` auto-detects protein (`.faa`) vs nucleotide (`.fna`, → Prodigal) input,
builds the databases on first run, then dispatches Snakemake.

## How it works

1. **Gene calls** — proteomes used directly; nucleotide assemblies → Prodigal.
2. **HMM scan** — `hmmscan` against `resources/hmm/ncycle_targets.hmm`, a
   concatenation of the **KOfam profile HMM for every KO in `config/targets.yaml`**
   + the **custom clade HMMs** (`amoA`, `nxrA`, `nxrB`, plus the octaheme HAO-family
   markers `hao` and `hdh` — trained CD-HIT→MAFFT→hmmbuild with calibrated trusted cutoffs) + a Pfam fallback for KO-less targets (archaeal
   amoA / PF12942). Thresholds (KOfam `ko_list`, custom TCs, per-KO overrides) live
   in `resources/hmm/tc_cutoffs.tsv`.
3. **BLAST gating** — `diamond blastp` against **curated, clade-spread UniProt seeds**
   (`resources/blast_db/`), per-target `blast_identity_min`, tagged `>{target_id}||{acc}`.
4. **Calls** — `apply_rules.py` integrates evidence per target. Signature
   precedence **custom HMM > KO > Pfam**; targets with
   `requires_blast_for_confirmation` (the homology traps) need a clade-specific
   BLAST hit or are `disqualified`.
5. **Reports** — per-sample `calls/ncycle_calls.tsv`, `complex_completeness.tsv`,
   `synergy_completeness.tsv`, `report/gap_analysis.txt`; cross-sample
   `multisample_matrix.tsv` + heatmap + per-pathway/complex/synergy/process figures.

## Pathways & target proteins

The tool resolves a genome's role across **all 9 nitrogen-cycle pathways** (49 marker
genes). Each gene is anchored on a KEGG Orthology (KO) where one exists; the KO column
below is what the KOfam backbone scans for (custom HMMs / Pfam noted where used).

### 1. Nitrogen fixation — N₂ → NH₃
| gene | KO | role |
|---|---|---|
| nifH | K02588 | nitrogenase Fe protein (dinitrogenase reductase) |
| nifD | K02586 | nitrogenase MoFe protein α |
| nifK | K02591 | nitrogenase MoFe protein β |
| nifE / nifN | K02587 / K02592 | FeMo-cofactor scaffold |
| nifB | K02585 | FeMo-cofactor biosynthesis (radical-SAM) |
| vnfH / vnfD | K22899 / K22896 | **V**-nitrogenase (alternative) Fe protein / α |
| anfG | K00531 | **Fe-only** nitrogenase δ (alternative) |

### 2. Nitrification — ammonia oxidation (NH₃ → NH₂OH → NO₂⁻)
| gene | KO | role |
|---|---|---|
| amoA | K28504 (custom HMM) | ammonia monooxygenase α — AOB / comammox |
| amoA_archaeal | PF12942 (Pfam-only) | ammonia monooxygenase α — **archaeal (AOA)** |
| amoB / amoC | K10945 / K10946 | ammonia monooxygenase β / γ |
| hao | K10535 (custom HMM) | hydroxylamine oxidoreductase |

### 3. Nitrification — nitrite oxidation (NO₂⁻ → NO₃⁻)
| gene | KO | role |
|---|---|---|
| nxrA | K00370 (custom HMM) | nitrite oxidoreductase α (catalytic, Mo) — NOB / comammox |
| nxrB | K00371 (custom HMM) | nitrite oxidoreductase β (Fe-S) |

### 4. Denitrification — NO₃⁻ → NO₂⁻ → NO → N₂O → N₂
| gene | KO | role |
|---|---|---|
| narG / narH / narI | K00370 / K00371 / K00374 | membrane nitrate reductase α / β / γ (cyt b) |
| napA / napB | K02567 / K02568 | periplasmic nitrate reductase catalytic / di-heme |
| nirK | K00368 | Cu nitrite reductase (NO-forming) |
| nirS | K15864 | cytochrome cd₁ nitrite reductase (NO-forming) |
| norB / norC | K04561 / K02305 | NO reductase large / small (cNOR) |
| norZ | K04748 | quinol-dependent NO reductase (qNOR) |
| nosZ | K00376 | nitrous oxide reductase (N₂O → N₂ sink) |

### 5. DNRA — dissimilatory nitrate/nitrite reduction to ammonium (retains N)
| gene | KO | role |
|---|---|---|
| nrfA | K03385 | cytochrome-c nitrite reductase, ammonia-forming (CCNiR) |
| nrfH | K15876 | NrfA-associated tetraheme cytochrome c |
| nirB / nirD | K00362 / K00363 | NADH nitrite reductase large / small |

### 6. Anammox — anaerobic ammonium oxidation (NH₄⁺ + NO₂⁻ → N₂)
| gene | KO | role |
|---|---|---|
| hzsA / hzsB / hzsC | K20932 / K20933 / K20934 | hydrazine synthase α / β / γ |
| hdh | K20935 (custom HMM) | hydrazine dehydrogenase (N₂H₄ → N₂) |

### 7. Assimilatory nitrate/nitrite reduction — NO₃⁻ → NH₄⁺ for biosynthesis
| gene | KO | role |
|---|---|---|
| narB | K00367 | ferredoxin-nitrate reductase |
| nasA / nasB | K00372 / K00360 | assimilatory nitrate reductase catalytic / electron-transfer |
| nirA | K00366 | ferredoxin-nitrite reductase |
| nasD | K17877 | assimilatory nitrite reductase (NAD(P)H) large |
| NR | K10534 | **eukaryotic** nitrate reductase (nia) |

### 8. Ammonia assimilation — NH₄⁺ → glutamate / glutamine
| gene | KO | role |
|---|---|---|
| glnA | K01915 | glutamine synthetase (GS) |
| gltB / gltD | K00265 / K00266 | glutamate synthase (GOGAT) large / small |
| gdhA | K00261/K00262 | glutamate dehydrogenase |

### 9. Organic-N mineralization — urea → NH₃
| gene | KO | role |
|---|---|---|
| ureC / ureB / ureA | K01428 / K01429 / K01430 | urease α (catalytic, Ni) / β / γ |
| ureG | K03190 | urease accessory (Ni insertion) |

### Obligatory complexes & process synergies

Beyond per-gene calls, the pipeline scores whether **multi-subunit complexes** are
complete (all subunits present) and whether **process modules** ("synergies") are
satisfied — this is what lets it say *"complete denitrifier"* vs *"has a stray nar gene"*.

- **Complexes (all subunits required):** nitrogenase (nifHDK), ammonia_monooxygenase
  (amoABC), nitrite_oxidoreductase (nxrAB), membrane_nitrate_reductase (narGHI),
  periplasmic_nitrate_reductase (napAB), NO_reductase_cNOR (norBC), hydrazine_synthase
  (hzsABC), glutamate_synthase (gltBD), urease (ureABC).
- **Synergies (process completeness):** complete_denitrification (narG+nirS+norB+nosZ),
  n2o_sink (nosZ), comammox (amoA+hao+nxrAB in one genome), anammox_complete
  (hzsABC+hdh), dnra_branch (napA+nrfA), nitrifier_denitrification (amoA+nirK+norB),
  assimilatory_complete (nasA+nirA+glnA+gltB).

## The homology traps — how they're resolved

`nxrA`/`narG` and `nxrB`/`narH` share KOfam KO **K00370/K00371**, so the KO tier
*cannot* separate them. The custom NOB-clade HMMs + BLAST gate do: `nxrA` scores
**F1 0.75 vs 0.18** for raw KOfam (which fires K00370 in every denitrifier), `nxrB`
**0.75 vs 0.22**. `amoA` (ammonia) is kept distinct from `pmoA` (methane) by an
AOB+comammox clade HMM: **F1 0.89 vs 0.73** (raw KOfam calls amoA on methanotroph
pmoA). Leave-one-genus-out CV confirms `nxrA` generalizes across a held-out
*Nitrospira* genus; `amoA` is clade-sensitive (comammox seeds included by design).

The one sequence-inseparable case — *Nitrobacter* NxrA ≈ NarG (~60% identity, scores
overlap) — is resolved by **operon synteny** on nucleotide/MAG input: a K00370/K00371 ORF
syntenic with narI (K00374, the respiratory-NAR cytochrome that NXR lacks) is narG/narH,
one without is nxrA/nxrB (`resolve_nxr_synteny` in `apply_rules.py`). Validated:
*N. winogradskyi* → both nxrA (standalone α) and narG (narGHJI operon); *E. coli* → narG
only, no nxr false-positive. Needs gene coordinates, so it is **nucleotide-input only**
(dormant for pre-called proteomes).

## Validated accuracy (reference panel: training + 8 hold-out genera)

The **de-leaked hold-out is the headline** — it is the only frame that estimates
generalization. The training panel is calibration data, so its perfect score is an
in-sample fit, not a generalization claim (read it as a rule-of-three per-cell error
bound, not "0 errors forever").

| Frame | precision | recall | micro-F1 | note |
|---|---|---|---|---|
| **Hold-out, de-leaked (HEADLINE)** | 0.86 [0.73, 0.96] | 0.82 [0.62, 0.97] | **0.84 [0.68, 0.95]** | 8 genera absent from training; seed-leaked cells removed |
| Hold-out, as-is (seed-contaminated) | 0.87 | 0.83 | 0.85 | before removing 6 same-strain/species seed cells |
| Whole panel (training + hold-out) | 0.97 [0.93, 0.99] | 0.96 [0.90, 1.00] | 0.96 [0.92, 0.99] | mixes in-sample + out-of-sample |
| Training panel (**in-sample fit**) | 1.00 | 1.00 | 1.00 | calibration data; error bound ≤0.47% (rule of 3) |

CIs are 95% genome-cluster percentile bootstrap (the genome, not the cell, is the unit of
independence). All 9 pathways have a validated positive. The de-leaked hold-out is the
credible generalization estimate; the training 1.00 must not be reported as such.

> **Why de-leaked?** A few hold-out cells were "detected" partly because a BLAST/HMM seed
> came from the hold-out organism's own strain/species (e.g. the γ-AOB seed D5BWX5 *is*
> *Nitrosococcus halophilus* Nc4, a hold-out genome). `validation/detect_seed_leakage.py`
> finds and removes these so the hold-out measures generalization, not memorization.

## Residuals — both now resolved

The two former residuals have been addressed; the only standing constraint is an
input-modality one (synteny needs gene coordinates):
- ***Nitrobacter* `nxrA/nxrB` ↔ `narG`** — NxrA is nested in the NAR clade (scores 543
  vs 544); no sequence method separates them. **Resolved by operon synteny** for
  nucleotide/MAG input (narI proximity; see above). On pre-called proteomes (no gene
  coordinates) it persists — use nucleotide/MAG input for NOB-vs-denitrifier resolution.
- **AOA `amoA_archaeal`** — was a gate misconfiguration, not a hard limit: the archaeal-
  specific PF12942 fires cleanly; the empty-seed BLAST gate was disqualifying it. Gate
  dropped → **F1 0.00 → 1.00**, ammonia-oxidation 0.93 → 0.98.

Everything else scores ≥0.94.

## Databases & reproducibility

`resources/.cache/profiles.tar.gz` is the one-time KOfam download (~1.5 GB; the
build extracts only the ~47 needed profiles). **The KOfam release is pinned by
content hash** (`KOFAM_PROFILES_SHA256` in `workflow/scripts/build_hmm_db.py`,
mirrored in `config/config.yaml`): the build verifies the tarball and writes
`resources/.cache/kofam_release.txt` as provenance. The validated baseline used
release **2026-05-24** (sha256 `b03d20b9…`); a mismatch warns loudly but does not
hard-fail. Deleting the cache reclaims space — it re-downloads on the next rebuild.

## Regression testing

`validation/test_regression.py` is an accuracy gate: it scores the panel and **fails
(exit 1) if any metric's 95% bootstrap CI lower bound drops below a locked-in floor**.
The floors are genome-cluster bootstrap CI lower bounds, not point estimates: training-panel
ALL micro-F1 and precision CI-lo ≥ 1.00 (the in-sample calibration scope, which scores
cleanly), homology-trap precision CI-lo ≥ 0.95, **de-leaked hold-out F1 CI-lo ≥ 0.65** (the
honest generalization guard), every pathway F1 CI-lo ≥ 0.70, plus point checks FP ≤ 4 and
scored cells ≥ 600. The generalization guarantee rests on the **de-leaked hold-out floor**,
not the in-sample 1.00. It shares its metric logic with `score_ncycle.py`
(`compute_metrics()`), so the gate and the report can never disagree.

```bash
make regression         # full: rebuild GT → run pipeline on $(PANEL) → score → gate
make regression-score   # fast: re-score existing results/ → gate (no pipeline run)
python validation/test_regression.py   # the gate alone (also runs under pytest)
```

A GitHub Actions workflow (`.github/workflows/regression.yml`) wires `make regression`
into CI with KOfam caching — it's a template pending two prerequisites: putting the repo
under git, and provisioning the reference panel (see the workflow comments).
