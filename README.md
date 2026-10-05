# ncycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **nitrogen cycle**.
Modeled on `Holomicrobiome-ewaste/ewaste-pipeline`; detection is **KO-primary**
(KOfam HMMs + adaptive per-KO thresholds) with custom clade HMMs and DIAMOND-BLAST
gating for the homology traps. Covers **51 targets** (47 KEGG Orthologies), **9
obligatory complexes**, and **11 process-completeness synergies** across all 9
N-cycle pathways.

**Status: validated.** End-to-end on the reference panel (training genomes + 8 independent
hold-out genera). The honest generalization headline is the **de-leaked hold-out micro-F1 ≈
0.84** (95% genome-cluster bootstrap CI [0.68, 0.95]; P 0.86, R 0.82). The **training-panel
F1 = 1.00 is an in-sample calibration fit** (thresholds/seeds tuned on those genomes; honest
per-cell error bound ≤0.47%, rule of three), **not** a generalization number. The
leakage-robust comparator advantage over raw KofamScan is **homology-trap precision**. See
[`validation/REPORT.md`](validation/REPORT.md) (§ Audit 2026-06-10 reframed the headline) and
[`ROADMAP.md`](ROADMAP.md).

## Install

Developed and validated on Linux (x86-64). Needs `git`, conda or mamba and `curl`;
about 3 GB of disk for the environment. Network access is needed for the install and
for the five genomes of the smoke test (NCBI); the databases are built from files in
the repository and the pipeline itself runs offline.

```bash
git clone https://github.com/Pedrolleao/ncycle-pipeline.git
cd ncycle-pipeline
conda env create -f envs/ncycle.yaml       # env `cycle-pipeline`; or: mamba env create …
conda activate cycle-pipeline
make test_protein     # fetch 5 genomes, build the databases, run end to end
make regression       # run the 39-genome reference panel (test_panel/) and check the accuracy floors
```

`make regression` ending in `OK: all … checks passed` means the install reproduces the
validated calls. `envs/ncycle.lock.yml` is the exact environment of the validation (linux-64), for
when the open version ranges of `envs/ncycle.yaml` resolve to something that behaves
differently.

## Run

```bash
python run.py --input <dir-of-.faa-or-.fna> --cores 8
# Outside the conda env, run.py re-runs itself inside `cycle-pipeline` (and creates it
# from envs/ncycle.yaml if it does not exist). To use another env with the same
# dependencies:
#   NCYCLE_ENV=<env-name> python run.py --input <dir>
```

`run.py` auto-detects protein (`.faa`) vs nucleotide (`.fna`, → Prodigal) input,
builds the databases on first run, then dispatches Snakemake. It asks whether
nucleotide input is isolate genomes or metagenome assemblies unless
`--prodigal-mode single|meta` is given, and writes the samples it found into the
`samples:` block of `config/config.yaml` — so `git status` shows that file as modified
after a run.

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
   `synergy_completeness.tsv`, `report/gap_analysis.txt`, `report/ncycle_map.*`;
   cross-sample `multisample_matrix.tsv`, figures, and an interactive `report.html`
   (see **Outputs** below).

## Outputs

All paths are under `paths.results_dir` (`results/` by default). Figures are written
as SVG (vector) and PNG (300 DPI).

**Per sample — `<sample>/`**

| file | what it is |
|---|---|
| `calls/ncycle_calls.tsv` | one row per target: status, evidence source, protein, HMM hit + E-value, BLAST reference + identity. For nucleotide input the columns `contig`, `start`, `end`, `strand` locate the called gene (empty for a pre-called proteome). One protein is reported per target; `other_copies` lists any further proteins that reach the same status (paralog copies) |
| `calls/complex_completeness.tsv`, `calls/synergy_completeness.tsv` | completeness of the 9 complexes and 11 process modules |
| `calls/ncycle_loci.tsv` | *(nucleotide input only)* called genes grouped into loci: genes on one contig with at most 5 other genes between them. `copy` says whether a gene is the one reported in `ncycle_calls.tsv` or an additional copy |
| `report/gap_analysis.txt` | plain-text summary |
| `report/ncycle_map.svg/.png` | the genome's calls drawn on the nitrogen cycle: each reaction arrow is solid (a complete route found), dashed (partial) or grey (absent), with the genes behind it |
| `report/loci.svg/.png` | *(nucleotide input only)* gene-arrow maps of every locus, grouped by pathway, to a common bp scale. Genes outside the target set are blank; a bar marks a contig end (where an operon may run off the assembly) |

**Across samples**

| file | what it is |
|---|---|
| `multisample_matrix.tsv` | genomes × (targets, complexes, modules) |
| `multisample_heatmap.svg/.png` | overview dot grid; genomes ordered by gene-content similarity |
| `figures/pathway_<pathway>.svg/.png` | one dot grid per pathway |
| `figures/complexes.svg/.png`, `figures/synergies.svg/.png` | complex / process-module completeness |
| `figures/ncycle_maps.svg/.png` | every genome's N-cycle map side by side (up to 48 genomes) |
| `report.html` | self-contained interactive report (no network needed): the gene grid and the complex / module grid with hover evidence, row search / ordering, and a per-genome panel with the N-cycle map, locus maps and the full calls table. Light and dark themes. Every figure in it (gene grid, complex / module grid, cycle map, each locus map) has a **Save PNG (300 dpi)** button: it downloads that figure as currently shown — row filter and order, hidden pathways, selected genome, light or dark theme — with its title and legend, rendered at 300 dpi (a grid too large for a browser canvas is saved at the highest resolution that fits, and says so). The page follows the group's *Simple Terminal* design system (`design/Simple`): JetBrains Mono, hairline `[ bracketed ]` frames, its dark palette or its Light variant according to the system theme, with a LIGHT / DARK selector in the top-right corner to pin either. The font is inlined from `workflow/scripts/fonts/` (SIL OFL 1.1, licence alongside), so the report looks the same offline and the PNG export uses it too; pathway colours stay the validated palette of the static figures. |

**Reading the glyphs** (same in every figure and in `report.html`): solid disc =
confirmed; half-filled = domain-only (HMM signature, no BLAST support); ring with a
cross = disqualified (failed the homology-trap gate); faint ring = absent. In the
complex / module grids: solid = complete, ring with `n/N` = partial, faint ring with
a cross = ruled out (an excluded gene is present, e.g. nosZ rules out `n2o_emitter`).
Colour always means pathway; the palette lives in `workflow/scripts/_domain.py`.

The figure and report scripts are shared, byte-identical, with the sulfur sister
pipeline; only `workflow/scripts/_domain.py` (palette, labels, file names) and
`workflow/scripts/_cycle_model.py` (the cycle diagram) are nitrogen-specific.

## Pathways & target proteins

The tool resolves a genome's role across **all 9 nitrogen-cycle pathways** (51 marker
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
| amoA_gamma | K28504 (custom HMM) | ammonia monooxygenase α — **γ-proteobacterial AOB** (*Nitrosococcus*) |
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
| nosZ | K00376 | nitrous oxide reductase, clade I (N₂O → N₂ sink) |
| nosZ_clade2 | K00376 (custom HMM) | nitrous oxide reductase, clade II / atypical (Sec-secreted) |

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
  assimilatory_complete (nasA+nirA+glnA+gltB), n2o_emitter (nir+nor present, no nosZ of
  either clade), n2o_sink_only (nosZ/nosZ_clade2 present, no nir/nor), nosZ_clade_I_likely
  (nosZ with nir/nor, no nosZ_clade2), nitrate_to_nitrite_leak (nar/nap present, no
  nirS/nirK/nrfA).

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

**The database inputs are pinned in the repository.** The 47 KOfam profiles and
their thresholds are built from `resources/kofam_pinned/` (release of 2026-05-24, the
one the tool was validated on), the Pfam fallback profiles from `resources/pfam_pinned/`
and the BLAST seeds from `resources/seeds_pinned/`. genome.jp serves KOfam from a
rolling URL and keeps no old releases — the release of 2026-09-29 has a different
threshold for 40 of the 47 KOs — and UniProt entries are revised and deleted, so
databases built from fresh downloads are not the validated ones. `build_hmm_db.py
--upstream` and `build_blast_db.py --upstream` download the current data anyway;
re-validate (`make regression`) before trusting the result.

**What a clone does not contain.** Pipeline results and the third-party tools and
databases of the comparator benchmark (METABOLIC, DRAM, NCycDB, GTDB proteomes). The scripts under `comparators/` and
the study configs `config/config_*.yaml` are the record of how the validation was run:
they carry paths of the machine it ran on and need editing to be re-run elsewhere.
Their outputs — the tables under `validation/` and `comparators/` — are committed.

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
make regression         # full: check the panel, rebuild GT → run pipeline on test_panel/ → score → gate
make regression-score   # fast: re-score existing results/ → gate (no pipeline run)
python validation/test_regression.py   # the gate alone (also runs under pytest)
```

A GitHub Actions workflow (`.github/workflows/regression.yml`) runs `make dbs`,
`make test_protein` and `make regression` on every push to `main`.
