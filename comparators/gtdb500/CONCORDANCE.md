# ncycle-pipeline — GTDB-scale Cross-Tool CONCORDANCE (pilot)

**Genomes:** 500 GTDB-representative (120 N-cycle-clade-enriched + 380 cross-phylum backbone) · **Tools:** ncycle, raw KofamScan, METABOLIC v4.0, NCycDB · **No ground truth** — agreement/divergence only.

> Concordance, not accuracy: arbitrary GTDB genomes have no curated truth. High agreement on specific markers + systematic divergence at the homology traps is the expected signature; the trap loci are where ncycle's clade resolution acts.


## 1. Pairwise agreement

Over targets both tools can represent. **Raw** = identical present/absent call (inflated by shared-absent cells). **Jaccard** = positive-call agreement |both present| / |either present| (the honest concordance metric). Reported for ALL targets and the homology-trap subset.

| tool pair | raw (ALL) | raw (trap) | Jaccard (ALL) | Jaccard (trap) |
|---|---|---|---|---|
| ncycle vs kofam | 0.828 | 0.804 | 0.290 | 0.218 |
| ncycle vs metabolic | 0.947 | 0.927 | 0.617 | 0.337 |
| ncycle vs ncycdb | 0.839 | 0.841 | 0.390 | 0.271 |
| kofam vs metabolic | 0.844 | 0.802 | 0.336 | 0.214 |
| kofam vs ncycdb | 0.779 | 0.758 | 0.358 | 0.303 |
| metabolic vs ncycdb | 0.838 | 0.813 | 0.405 | 0.224 |

## 2. nxrA↔narG trap locus (headline)

- ncycle calls **nxrA** in **18** genomes.
- Of those, **12** have ≥1 comparator calling **narG (and NOT nxrA)** on the same genome — the mis-routing of nitrite-oxidizer NxrA into the nitrate-reductase family that the curated-panel benchmark documented, now observed at scale.


## 3. Targets where ncycle is DISTINCTIVE (differs from every comparator)

| target | n genomes ncycle-distinct | trap? |
|---|---|---|
| nasD | 130 |  |
| norZ | 127 |  |
| gltD | 106 |  |
| norC | 54 |  |
| vnfH | 48 |  |
| narG | 40 | trap |
| norB | 38 |  |
| narH | 22 | trap |
| amoA_gamma | 21 |  |
| narB | 16 |  |
| nirD | 15 | trap |
| nxrB | 11 | trap |
| nosZ_clade2 | 11 |  |
| vnfD | 5 |  |
| nasA | 2 |  |

## 4. Per-target present-call counts (by tool; '-' = tool cannot represent)

Sorted by divergence (max−min present-count among representing tools).

| target | ncycle | kofam | metabolic | ncycdb | divergence |
|---|---|---|---|---|---|
| norZ | 26 | 497 | 153 | - | 471 |
| gdhA | 232 | 48 | 232 | 445 | 397 |
| nasA | 45 | 211 | 31 | 340 | 309 |
| gltD | 159 | 344 | 271 | 455 | 296 |
| nirB | 50 | 5 | 43 | 263 | 258 |
| nasB | 23 | 255 | 1 | 56 | 254 |
| napA | 38 | 273 | 33 | 150 | 240 |
| NR | 8 | 12 | 0 | 230 | 230 |
| nosZ | 40 | 72 | 38 | 260 | 222 |
| gltB | 144 | 3 | 143 | 210 | 207 |
| nxrA | 18 | 208 | 43 | 5 | 203 |
| norC | 37 | 220 | 95 | 114 | 183 |
| nirD | 13 | 33 | 51 | 194 | 181 |
| nirA | 48 | 100 | 22 | 199 | 177 |
| narG | 36 | 208 | 43 | 76 | 172 |
| nirS | 18 | 58 | 17 | 183 | 166 |
| nirK | 106 | 173 | 75 | 206 | 131 |
| nrfH | 132 | 178 | 47 | - | 131 |
| nasD | 130 | 0 | 0 | - | 130 |
| hao | 94 | 218 | 93 | 112 | 125 |
| ureC | 62 | 10 | 68 | 135 | 125 |
| narB | 53 | 125 | 1 | 106 | 124 |
| nxrB | 19 | 140 | 36 | 33 | 121 |
| nifB | 47 | 155 | 47 | - | 108 |
| narH | 39 | 140 | 36 | 99 | 104 |
| norB | 43 | 139 | 81 | 93 | 96 |
| hzsA | 3 | 77 | 5 | 21 | 74 |
| ureA | 64 | 2 | 66 | 70 | 68 |
| nifK | 29 | 47 | 26 | 92 | 66 |
| ureB | 61 | 3 | 64 | 66 | 63 |
| nifH | 40 | 93 | 46 | 100 | 60 |
| nrfA | 40 | 92 | 75 | 42 | 52 |
| vnfH | 50 | 4 | 0 | - | 50 |
| narI | 43 | 90 | 43 | 70 | 47 |
| nosZ_clade2 | 28 | 72 | 38 | - | 44 |
| nifD | 28 | 7 | 27 | 50 | 43 |
| amoB | 21 | 50 | 13 | 38 | 37 |
| amoC | 35 | 46 | 18 | 41 | 28 |
| amoA_gamma | 8 | 15 | 32 | - | 24 |
| glnA | 449 | 449 | 449 | 467 | 18 |
| amoA | 16 | 15 | 32 | 21 | 17 |
| nifN | 27 | 39 | 26 | - | 13 |
| hzsC | 5 | 14 | 5 | 7 | 9 |
| hdh | 8 | 0 | 8 | 9 | 9 |
| napB | 23 | 23 | 28 | 20 | 8 |
| hzsB | 7 | 9 | 3 | 11 | 8 |
| vnfD | 6 | 0 | 1 | - | 6 |
| anfG | 3 | 0 | 2 | 2 | 3 |
| nifE | 29 | 30 | 29 | - | 1 |
| amoA_archaeal | 20 | - | - | 20 | 0 |
| ureG | 67 | 67 | 67 | - | 0 |
