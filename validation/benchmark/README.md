# Benchmark harness — ncycle-pipeline vs external N-cycle tools

Statistically-robust head-to-head on the shared training panel (31 genomes, `curated_v12`
ground truth; hold-outs excluded as pre-registered). Design and decision rules are
pre-registered in [`prereg.md`](prereg.md) (+ 2026-06-09 amendment adding NCycDB).

## Files
- `prereg.md` — pre-registration (hypotheses, comparators, vocabulary, statistics, decisions).
- `adapters.py` — loads each tool's calls into one shared vocabulary + ground truth.
  Working: `ncycle` (our matrix), `kofam` (raw KO baseline). Scaffolded: `metabolic`,
  `ncycdb`, `dram` (parse a normalized TSV → our targets; return nothing until provided).
- `benchmark_stats.py` — the statistics engine (genome cluster bootstrap CIs, paired Δ
  bootstrap, exact McNemar, BH-FDR; subunit + step resolutions).
- `benchmark_results.tsv` — last run's full table.

## Run
```bash
# whatever is available (ncycle + kofam need no extra input):
python validation/benchmark/benchmark_stats.py --bootstrap 10000

# add an external tool once you have its normalized output:
python validation/benchmark/benchmark_stats.py \
    --metabolic metabolic.tsv --ncycdb ncycdb.tsv --dram dram.tsv
```

## Statistics (why it's robust)
- **Unit of independence = genome** (cells within a genome are correlated) — all 95% CIs
  use a **genome-level cluster bootstrap** (B=10,000), not naive per-cell resampling.
- **Paired Δ bootstrap** (ours − comparator on the same resampled genomes) is the primary
  significance test; a win requires the Δ CI to exclude 0.
- **Exact McNemar** on cell-level discordant pairs is reported as a secondary check.
- **Benjamini–Hochberg FDR** corrects the pre-specified family {comparator}×{trap-precision, ALL-F1}.
- Two resolutions reported: **subunit** (49 targets; the homology-trap claim) and **step**
  (13 transformations; fairest to coarse tools). The `compute_metrics()` math is unit-tested.

## Adding a comparator (3 steps)
1. Run the tool on the panel with **default settings**; record version + command in `adapters.py`.
2. Export its native output to a normalized TSV (see each loader's docstring for columns)
   and complete the name→target map (`NCYC_MAP` / `DRAM_STEP_MAP`, or the KO export for METABOLIC).
3. Re-run `benchmark_stats.py --<tool> PATH`. The new tool joins the CI/Δ/McNemar/BH tables.

## Head-to-head results (curated_v12 panel, 31 training genomes, B=10,000)
**All four comparators run: raw KofamScan, METABOLIC v4.0, NCycDB, DRAM v1.4.6.** Both
pre-registered confirmatory endpoints survive BH correction against *every* comparator (8/8) →
**ncycle WINS**. Full narrative + per-target confusions in [`COMPARISON_REPORT.md`](COMPARISON_REPORT.md).

| subunit metric | ncycle | KofamScan | METABOLIC | NCycDB | DRAM |
|---|---|---|---|---|---|
| ALL micro-F1 | **1.000** | 0.916 | 0.887 | 0.858 | 0.544 |
| **homology-trap precision** | **1.000** | 0.836 | 0.810 | 0.791 | 0.437 |
| homology-trap recall | **1.000** | 0.850 | 0.567 | 0.883 | 0.983 |

Each comparator fails the homology traps differently: **KofamScan over-calls** (shared KOs → all
targets), **METABOLIC under-calls** (stricter KOfam thresholds, trap recall 0.57), **NCycDB
mis-routes** (curated nxrA family exists but DIAMOND best-hit sends NOB nxrA → narG; nitrite-ox
recall 0), **DRAM coarsely over-calls** (module expansion → trap precision 0.44; also no module for
glnA/ureC). ncycle's clade-HMM + BLAST-gate + synteny layers resolve all three trap classes →
trap precision 1.00. NCycDB is the strongest comparator (alone resolves archaeal amoA, avoids the
methanotroph FP); DRAM the weakest (module-level can't address subunit traps).
