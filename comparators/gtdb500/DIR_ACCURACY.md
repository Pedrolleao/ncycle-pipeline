# ncycle — nxrA-vs-narG identity accuracy vs phylogeny-anchored reference (B9)

Orthogonal, **non-circular** check: ncycle disambiguates the shared-KO nxrA/narG trap with custom HMMs + gating; the reference call comes from type-II DMSO-reductase sequence ancestry (type-labeled nxr/nar refs, best-hit DIAMOND). Agreement = accuracy of the identity call. Unit = each nxrA/narG call.

**nxrA/narG calls scored:** 54 · **with a reference hit:** 52

| stratum | agreement | accuracy | Wilson 95% CI |
|---|---|---|---|
| **Characterized (named genus)** | 20/28 | **0.714** | [0.529, 0.847] |
| Uncultured candidate phyla (`p__`) | 17/24 | 0.708 | [0.508, 0.851] |
| Confident best-hit (margin ≥ 0.10) | 26/37 | 0.703 | [0.542, 0.825] |
| **All with a hit** | 37/52 | 0.712 | [0.577, 0.817] |

> **Interpretation — best-hit is INCONCLUSIVE for the nxr/nar trap (this is the finding, not an ncycle accuracy verdict).** Unlike the dsrAB reductive/oxidative split (cleanly separable, scycle validated at 100% on characterized genomes by the same method), NxrA and NarG are reciprocally close paralogs in the type-II DMSO-reductase family: Nitrobacter-type NXR is phylogenetically embedded next to NarG. The discordances concentrate on exactly the expected lineages — **Nitrobacter ×5, Nitrococcus, Nitrolancea** (NOB) calling narG but best-hitting the nearby Nitrobacter NxrA reference (no Nitrobacter-clade narG exists in the 20-seq reference), plus divergent candidate phyla. In these, ncycle's custom-HMM+gating call is the more reliable signal and best-hit is mislabeling. The nxr/nar identity therefore REQUIRES the rigorous placement tier (B9.4, with operon/synteny context); the curated-panel benchmark remains the N accuracy authority. This asymmetry (sulfur resolvable by best-hit, nitrogen not) is itself a reportable result.


## Confusion matrix — ncycle call (rows) × phylogeny best-hit (cols)

| ncycle ↓ \ phylo → | nxr | nar | no_hit |
|---|---|---|---|
| **nxrA** | 12 | 6 | 0 |
| **narG** | 9 | 25 | 2 |

## Escalation set for B9.4 (EPA-ng/gappa placement)

- **Disagreements** (ncycle ≠ phylo, with hit): **15**.
- **Low best-hit margin** (<0.10): **15** — Nitrobacter-type NxrA / divergent periplasmic NXR sit close to NarG; place these.
- **No reference hit**: **2** (divergent enzyme → placement needed).

Per-call detail: `DIR_ACCURACY.tsv`. Lightweight tier; rigorous placement is B9.4.

