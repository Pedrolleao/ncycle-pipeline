# GTDB-scale cross-tool benchmark — 500-genome report

**Status:** complete (500 GTDB r232 representative genomes). This is the scaled-up run that the
111-genome pilot (`../gtdb_pilot/PILOT_REPORT.md`) was sized to budget. Same locked design:
a **concordance** study — arbitrary GTDB genomes have no curated ground truth, so we measure
where tools **agree/disagree**, not accuracy. Comparators: **ncycle, raw KofamScan,
METABOLIC v4.0, NCycDB** (DRAM excluded at scale, per prereg). The accuracy verdict remains the
645-cell curated-panel benchmark (`../../validation/benchmark/COMPARISON_REPORT.md`); this run
shows that the trap divergence it documented is **pervasive across GTDB**.

---

## 1. Panel
500 GTDB r232 species-representatives (`select_gtdb_pilot.py`, seed=1234, N_BACKBONE raised):
- **120 N-cycle-clade-enriched** across the trap-relevant clades (NOB *Nitrospira/Nitrobacter/
  Nitrococcus/Nitrospina/Nitrolancea*; AOB *Nitrosomonas/Nitrosospira/Nitrosococcus*; AOA
  *Nitrosopumilus/Nitrososphaera/Nitrosarchaeum*; anammox Brocadiales; methanotroph
  Methylococcaceae/Methylacidiphilaceae; denitrifiers).
- **380 cross-phylum backbone** (random species-reps spanning broad GTDB phylogenetic breadth).
- Fetched via NCBI `datasets`, Prodigal-called uniformly — all four tools see identical proteins.

Trap loci are well-populated by ncycle present-calls: nxrA 18, narG 36, amoA 16,
amoA_archaeal 20, amoA_gamma 8, nosZ 40, nirK 106, hzsA 3.

---

## 2. Wall-clock — measured on this 32-core box

| stage | 500 genomes | notes |
|---|---|---|
| ncycle full pipeline (incl. hmmscan over the KO HMM DB) | **~27 min** | proteomes→matrix; not the bottleneck |
| NCycDB (DIAMOND vs NCyc_100) | **~25 min** | |
| METABOLIC v4.0 (`-kofam-db full`) | **hours — the binding constraint** | KOfam hmmsearch ~3 h + KEGG module/identifier downstream |

**Scaling finding holds at 500:** ncycle's hmmscan parallelizes cleanly (~27 min for 500), so it
is *not* the GTDB-scale bottleneck. The cost is entirely **METABOLIC's full-KOfam hmmsearch +
its downstream KEGG-module/identifier passes**, consistent with the pilot's ~6 h/500 budget.

**Operational note (METABOLIC resume):** the METABOLIC run was interrupted partway and resumed
without recomputing the 3 h hmmsearch, via a shim (`_shim/hmmsearch`) that short-circuits any
search whose `--tblout` already exists and is non-empty, re-running only the uncached steps. The
run was stopped at the **KEGG-identifier step** — the only output the benchmark needs — skipping
dbCAN2/MEROPS/figure generation. All **500 `KEGG_identifier_result/*.result.txt`** were verified
present before normalization. (`_shim/kegg.status` shows a `TIMEOUT` from the auto-watcher, which
gave up before the KEGG step began; the perl run nevertheless completed that step at 11:47:50,
which is what matters.)

---

## 3. Concordance findings (full tables in `CONCORDANCE.md`)

### 3.1 Pairwise agreement
Raw present/absent agreement is high but **inflated by shared-absent cells**; the **Jaccard
positive-call agreement** is the honest metric.

| tool pair | Jaccard (ALL) | Jaccard (**trap**) |
|---|---|---|
| ncycle vs KofamScan | 0.290 | 0.218 |
| ncycle vs METABOLIC | 0.617 | 0.337 |
| ncycle vs NCycDB | 0.390 | 0.271 |
| KofamScan vs METABOLIC | 0.336 | 0.214 |
| KofamScan vs NCycDB | 0.358 | 0.303 |
| METABOLIC vs NCycDB | 0.405 | 0.224 |

**Headline:** at the homology-trap loci **no two tools agree on even 60 % of positive calls; most
pairs sit at 0.21–0.34.** Tool choice materially changes which N-cycle genes get reported exactly
where the biology is hardest. The two **threshold-gated** tools — ncycle and METABOLIC (METABOLIC
applies KOfam adaptive thresholds) — are the most similar **overall** (Jaccard 0.617) yet still
diverge sharply at the **traps** (0.337). **Raw KofamScan over-calls broadly** and is the divergent
KO tool (ncycle vs KofamScan 0.290; KofamScan vs METABOLIC 0.336): with no thresholds it fires every
borderline KO hit (nxrB 52 vs METABOLIC 36; amoB 50 vs 13; amoC 47 vs 18). NCycDB is the other
global outlier (Jaccard ~0.36–0.41 with everyone) — loose DIAMOND best-hit over-calls broadly
(nasA 340/500, nirB 263, gdhA 445).

> **Numbers corrected 2026-06-15:** an earlier version of this table (pre-B8 raw-KofamScan loader)
> showed ncycle vs KofamScan 0.673 / KofamScan vs METABOLIC 0.798. The B8 fix to
> `adapters.load_kofam` (KO/score columns for the hmmsearch domtblout layout) was applied but the
> GTDB concordance was not re-run until now. The fix only moves the **raw-KofamScan** pairs (every
> non-kofam pair is byte-identical), and it sharpens the story: raw KofamScan over-calls, so the two
> threshold-gated tools (ncycle, METABOLIC) are now correctly the most concordant.

### 3.2 nxrA↔narG mis-routing reproduced at scale
ncycle calls **nxrA** in **18** genomes; in **12** of them ≥1 comparator calls **narG (not nxrA)**
on the same genome — the nitrite-oxidizer NxrA→NarG mis-routing the curated-panel benchmark
documented, now observed across 500 GTDB representatives.

### 3.3 ncycle-distinctive loci (differs from every representing comparator)
Largest: **nasD 130, norZ 127, gltD 126, nrfH 74, norC 58, amoA_gamma 36** — dominated by
BLAST-gated targets the KO-only tools cannot separate (nasD/nrfH have no clean KO; norZ vs norB is
the shared-K04561 over-call by the KO tools) plus the nxr/nar/amo trap subunits. The pattern
matches the curated-panel mechanisms exactly.

### 3.4 Comparator coverage gaps (documented, not divergence)
METABOLIC's bundled KOfam lacks K10534 (NR) and K17877 (nasD) → both 0/500 for METABOLIC (≤2-cell
known gap). NCycDB cannot represent BLAST-only targets (nasD, nrfH, norZ, vnfH, amoA_gamma,
nosZ_clade2, several nif accessories) — reported as `-`, not as disagreement.

---

## 4. Honest caveats
- **Concordance ≠ accuracy.** Without ground truth this measures *divergence*, not who is correct.
  For the trap loci the curated-panel benchmark (645-cell, bootstrap-CI, BH-corrected) already
  established ncycle as correct; this run shows the divergence is pervasive at scale.
- For widespread housekeeping-ish targets (gdhA, gltD, glnA) the higher NCycDB counts may be
  biologically correct rather than over-calls; we do not adjudicate those here.
- Enrichment means the panel is **not** a random GTDB draw (by design, to populate trap loci);
  agreement numbers are per-genome call concordance, not population frequencies.

---

## 5. Reproduction
```bash
# selection (seed=1234, N_BACKBONE raised to reach 500) → selection.tsv / accs.txt / panel_map.tsv
python3 comparators/gtdb_pilot/select_gtdb_pilot.py
# datasets download → gtdb_dl/ ; prodigal → proteomes/
snakemake --cores 30 --configfile config/config_gtdb_pilot.yaml --snakefile workflow/Snakefile  # ncycle + kofam
python3 comparators/build_ncycdb_tsv.py   --proteome-dir comparators/gtdb500/proteomes \
        --out comparators/gtdb500/ncycdb.tsv --outdir comparators/gtdb500/ncycdb_out
# METABOLIC-G.pl -t 32 -in metabolic_in -kofam-db full -o metabolic_out   (resume via _shim/ if interrupted)
python3 comparators/build_metabolic_tsv.py --kegg-dir comparators/gtdb500/metabolic_out/KEGG_identifier_result \
        --out comparators/gtdb500/metabolic.tsv
python3 validation/benchmark/concordance.py --results comparators/gtdb500/results \
        --ncycdb comparators/gtdb500/ncycdb.tsv --metabolic comparators/gtdb500/metabolic.tsv \
        --panel-map comparators/gtdb500/panel_map.tsv --out comparators/gtdb500/CONCORDANCE.md
```

---

## 6. Outcome
The harness ran end-to-end at 500 genomes and the concordance signal is clear, on-message, and a
faithful scale-up of the 111-genome pilot: **pervasive trap divergence (no pair >0.6 Jaccard on
trap positives), the nxrA→narG mis-routing reproduced (12/18), and ncycle distinctive precisely on
the BLAST-gated and trap targets.** METABOLIC remains the sole scaling constraint (~hours);
dropping it would let ncycle + NCycDB + KofamScan reach 1,000 genomes in <1 h if breadth beyond 500
is wanted. No further work is required to close this benchmark line.
