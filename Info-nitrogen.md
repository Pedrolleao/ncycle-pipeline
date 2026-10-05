# Marker Genes & Proteins of the Microbial Nitrogen Cycle

> A pathway-by-pathway atlas of the nitrogen cycle and every protein involved in each step, with the marker identifiers needed to detect them in MAGs. Built as the reference layer for a MAG-mapping tool modeled on the **Holomicrobiome-ewaste** pipeline (see `targets.yaml` for the tool-ready inventory derived from this document).
>
> **Identifier backbone.** Each enzyme is anchored on **KEGG Orthology (KO)** numbers first — the nitrogen cycle is exhaustively curated in KEGG — then on **Pfam / InterPro** and **TIGRFAM / NCBIfam** HMM profiles, with canonical reference organisms for UniProt seed curation.
>
> **⚠️ convention.** Entries marked ⚠️ share their KO and/or Pfam with a different reaction (homology traps). For these, KO/Pfam alone over-calls — they require BLAST confirmation against curated seeds (the "dsrA-style gate" used throughout the ewaste pipeline) and, in the tool, a calibrated custom HMM. The three load-bearing traps in the N-cycle are: **`nxrAB` ↔ `narGH`** (identical KOs K00370/K00371), **`amoA` ↔ `pmoA`** (ammonia vs methane monooxygenase, shared B/C subunits), and **`nirB`/`nirA`/`nirD` ↔ `dsrAB`** (shared PF01077 NIR_SIR fold).
>
> Always cross-check at [KEGG](https://www.genome.jp/kegg/), [InterPro/Pfam](https://www.ebi.ac.uk/interpro/), and [NCBI NCBIfam](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/) before building HMM searches.

---

## The nitrogen cycle at a glance

```
                      N2 (dinitrogen gas)
                   ┌────────┘   ↑      ↑
       N2 fixation │            │      │ denitrification
        (nif/vnf)  ↓            │      │ terminal step (nosZ)
                  NH3/NH4+      │     N2O
        ┌──────────┤           │      ↑ (nor: norBC)
 assimilation      │ nitrification    NO
 (GS/GOGAT)        │ step 1 (amo,hao)  ↑ (nir: nirK/nirS)
        │          ↓                  NO2-
   organic-N    NH2OH → NO2- ←────────┤
   (urease) ←──┘        │             │ denitrification / DNRA
        ↑               │ nitrification│ (nar/nap → nir/nrf)
   mineralization       │ step 2 (nxr) ↓
                        └────────→ NO3-
                              ↑    │
              assimilatory ───┘    └─── DNRA → NH4+ (nrfA)
              (nas/narB/nirA)            anammox: NH4+ + NO2- → N2 (hzs, hdh)
```

Eight transformations move nitrogen between five oxidation states (NH₃ −3 → N₂ 0 → N₂O +1 → NO +2 → NO₂⁻ +3 → NO₃⁻ +5). Sections 1–10 below document each, then the obligatory-complex and process-completeness sections interpret co-occurrence.

---

## 1. Nitrogen Fixation (N₂ → NH₃)

Reduction of atmospheric N₂ to bioavailable ammonia by the nitrogenase complex — the sole biological entry point of new fixed nitrogen. KEGG module **M00175**.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Nitrogenase iron protein (dinitrogenase reductase) | *nifH* | K02588 | PF00142 | IPR005977 | TIGR01287 | Fer4_NifH (ATP-binding Fe-protein) | N₂ fixation — electron delivery | The universal nitrogenase marker and phylogenetic backbone for diazotroph surveys; delivers ATP-driven electrons to the MoFe protein. Presence of *nifH* is the single best proxy for fixation capability |
| Nitrogenase MoFe protein α subunit | *nifD* | K02586 | PF00148 | IPR000510 | TIGR01282 | Oxidored_nitro (MoFe α) | N₂ fixation — catalysis | Houses the FeMo-cofactor where N₂ is reduced; obligate partner of NifK. *nifD/nifK* discriminate true Mo-nitrogenase from alternative systems |
| Nitrogenase MoFe protein β subunit | *nifK* | K02591 | PF00148 | IPR000510 | TIGR01286 | Oxidored_nitro (MoFe β) | N₂ fixation — catalysis | Forms the α₂β₂ catalytic tetramer with NifD; non-functional alone |
| Nitrogenase cofactor scaffold / maturation | *nifE, nifN, nifB* | K02587, K02592, K02585 | PF00148 / PF06180 | — | TIGR01283 (nifE) | FeMo-co biosynthesis | N₂ fixation — cofactor assembly | Required to build FeMo-co; their presence distinguishes a complete, functional fixation system from orphan *nifH* hits |
| V-nitrogenase (alternative) | *vnfD/G/K, vnfH* | K22896, K22897, K22898, K22899 | PF00148 / PF00142 | — | — | Vanadium nitrogenase | N₂ fixation — Mo-independent | Used under Mo limitation; flags diazotrophs with alternative-nitrogenase repertoire |
| Fe-only nitrogenase δ / alt. | *anfG (δ)* | K00531 | PF03139 | — | — | Nitrogenase δ subunit | N₂ fixation — Mo/V-independent | Third nitrogenase isoform; least efficient, niche-specific |

---

## 2. Nitrification — Step 1: Ammonia Oxidation (NH₃ → NH₂OH → NO₂⁻)

Aerobic oxidation of ammonia, performed by ammonia-oxidizing bacteria (AOB), archaea (AOA), and the first half of comammox. KEGG module **M00528**.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Ammonia monooxygenase subunit A | *amoA* ⚠️ | K28504 (amoA KO); **K10944** is the *pmoA* (methane) KO | PF02461 (bacterial); **PF12942** (archaeal AmoA) | IPR003393 | TIGR03080 (pmoA/amoA) | Copper membrane monooxygenase, subunit A | NH₃ → NH₂OH (catalytic) | The defining nitrification marker and rate-limiting first step. ⚠️ **Shares family with methane monooxygenase `pmoA`** — and K28504 itself is NOT a clean ammonia/methane discriminator: methanotroph *pmoA* (both verrucomicrobial and γ-proteobacterial) and archaeal *amoA* also fire K28504 at scores near its KOfam threshold. Distinguishing AOB *amoA*, AOA *amoA*, comammox *amoA*, γ-AOB *amoA*, and methanotroph *pmoA* requires phylogeny/BLAST gating + clade HMMs, not the KO |
| Ammonia monooxygenase subunit B | *amoB* ⚠️ | K10945 (*pmoB-amoB*) | PF04744 | IPR006833 | — | Monooxygenase_B | NH₃ → NH₂OH (structural) | ⚠️ Subunit B is shared verbatim with methanotroph *pmoB* at the KO level |
| Ammonia monooxygenase subunit C | *amoC* ⚠️ | K10946 (*pmoC-amoC*) | PF04896 ⚠️ | IPR006980 | — | PmoC/AmoC | NH₃ → NH₂OH (structural) | ⚠️ Shared with *pmoC*; completes the AmoCAB membrane complex |
| Hydroxylamine oxidoreductase / dehydrogenase | *hao* ⚠️ | K10535 | ⚠️ multiheme cyt c clan (no dedicated Pfam) | IPR036866 (rel.) | — | Octaheme cytochrome c (P460) | NH₂OH → NO/NO₂⁻ | Oxidizes hydroxylamine, recovering electrons. ⚠️ Octaheme fold is shared with NrfA-family multiheme cytochromes — **KO-primary + custom HMM/BLAST detection** recommended (no clean Pfam) |

> **AOA note.** Archaeal ammonia oxidizers (*Nitrosopumilus*, *Nitrososphaera*) use a divergent AmoA (Pfam **PF12942**) and lack a canonical *hao*; their hydroxylamine handling differs. Treat archaeal *amoA* as a distinct target from bacterial *amoA*.

---

## 3. Nitrification — Step 2: Nitrite Oxidation (NO₂⁻ → NO₃⁻)

Oxidation of nitrite to nitrate by nitrite-oxidizing bacteria (NOB: *Nitrobacter*, *Nitrospira*, *Nitrococcus*) via nitrite oxidoreductase (NXR).

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Nitrite oxidoreductase α (catalytic) | *nxrA* ⚠️ | **K00370** (identical to *narG*) | PF00384 + PF01568 | IPR006468 | TIGR01580 | Molybdopterin oxidoreductase | NO₂⁻ → NO₃⁻ | Completes nitrification; NXR runs the same chemistry as nitrate reductase in reverse. ⚠️ **`nxrA` and `narG` share KO K00370 and Pfam architecture** — the single most important N-cycle disambiguation. Requires BLAST gating against NOB NxrA seeds + operon/phylogeny context |
| Nitrite oxidoreductase β (Fe-S) | *nxrB* ⚠️ | **K00371** (identical to *narH*) | PF13247 / PF13534 (4Fe-4S) | IPR017896 | TIGR01660 | Beta/Fe-S subunit | NO₂⁻ → NO₃⁻ | ⚠️ Shares KO K00371 with *narH*; *nxrB* is often used as the preferred NOB phylomarker because it is somewhat more clade-resolved than *nxrA* |

---

## 4. Complete Nitrification — Comammox (NH₃ → NO₂⁻ → NO₃⁻ in one organism)

*Nitrospira* lineage II clades that carry **both** ammonia-oxidation and nitrite-oxidation machinery. KEGG module **M00804**. No new genes — comammox is detected as **co-occurrence** of `amoCAB` + `hao` + `nxrAB` in a single MAG, with a distinct comammox-clade *amoA*.

| Component | Genes | KO set | Relevance |
|---|---|---|---|
| Ammonia oxidation | *amoCAB*, *hao* | K28504, K10945, K10946, K10535 | Comammox *amoA* forms its own phylogenetic cluster, separable from canonical AOB/AOA *amoA* by sequence — the key marker for comammox identification |
| Nitrite oxidation | *nxrAB* | K00370, K00371 | Same NXR as section 3; in comammox it belongs to the *Nitrospira* NXR clade |

> Comammox completeness is a **module-level call** (see Process-Completeness section), not a single gene.

---

## 5. Denitrification (NO₃⁻ → NO₂⁻ → NO → N₂O → N₂)

Stepwise anaerobic respiration of nitrogen oxides, returning fixed N to the atmosphere. Often modular across taxa (partial denitrifiers). KEGG module **M00529**.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Membrane nitrate reductase α | *narG* ⚠️ | K00370 | PF00384 + PF01568 | IPR006468 | TIGR01580 | Molybdopterin oxidoreductase | NO₃⁻ → NO₂⁻ | Dissimilatory nitrate reduction (first denitrification step). ⚠️ Shares KO/Pfam with *nxrA* — BLAST-gate to denitrifier NarG seeds |
| Membrane nitrate reductase β | *narH* ⚠️ | K00371 | PF13247 / PF13534 | IPR017896 | TIGR01660 | Fe-S subunit | NO₃⁻ → NO₂⁻ | ⚠️ Shares K00371 with *nxrB* |
| Membrane nitrate reductase γ | *narI* | K00374 | PF02665 | IPR003816 | TIGR00351 | Nitrate reductase γ (cyt b) | NO₃⁻ → NO₂⁻ | Membrane anchor; completes the NarGHI complex |
| Periplasmic nitrate reductase | *napA* | K02567 | PF00384 + PF01568 | IPR006657 | TIGR00509 | Molybdopterin oxidoreductase | NO₃⁻ → NO₂⁻ | Periplasmic alternative to Nar; functions aerobically/microaerobically. Also shared with DNRA. ⚠️ Mo-oxidoreductase family is broad — BLAST-gate (see ewaste `napA` precedent) |
| Periplasmic nitrate reductase, small | *napB* | K02568 | PF03892 | IPR005591 | — | NapB di-heme cytochrome c | NO₃⁻ → NO₂⁻ | Electron-donor partner of NapA |
| Cu-type nitrite reductase | *nirK* | K00368 | PF00394 + PF07731 + PF07732 | IPR001287 | TIGR02376 | Multicopper oxidase (Cu-NIR) | **NO₂⁻ → NO** (committing step) | Produces nitric oxide — the committed, defining step of denitrification (the gas-producing branch point vs DNRA). Copper-dependent; mutually (largely) exclusive with *nirS* in a genome |
| cd₁-type nitrite reductase | *nirS* | K15864 | PF13442 ⚠️ + PF02239 | IPR002572 | TIGR02375 | Cytochrome cd₁ | **NO₂⁻ → NO** (committing step) | Heme cd₁ alternative to NirK; same reaction, distinct fold |
| Nitric oxide reductase, large (cNOR) | *norB* | K04561 | PF00115 | IPR000883 | — | Heme-copper oxidase (NO red.) | NO → N₂O | Detoxifies/respires NO; large subunit of cNOR. qNOR (*norZ*, single-subunit) is an alternative |
| Nitric oxide reductase, small (cNOR) | *norC* | K02305 | PF00034 | IPR009056 | — | Cytochrome c (NorC) | NO → N₂O | Electron-donor subunit of cNOR |
| Nitrous oxide reductase | *nosZ* ⚠️ | K00376 | PF18764 (catalytic) + PF00116 (CuA) ⚠️ | IPR026434 | TIGR04244 | Multicopper (CuZ/CuA) N₂O reductase | **N₂O → N₂** (terminal step) | The only biological N₂O sink; its presence/absence determines whether a denitrifier is a net N₂O **source or sink** (greenhouse-gas relevance). Clade I and Clade II (*nosZ*-II often without *nir/nor*) are ecologically distinct |

---

## 6. Dissimilatory Nitrate Reduction to Ammonium — DNRA (NO₃⁻ → NO₂⁻ → NH₄⁺)

Respiratory route that **retains** nitrogen as ammonium instead of losing it as gas — the competitor of denitrification at the nitrite branch point. KEGG module **M00530**.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Nitrate reductase (shared) | *narG/H/I* or *napA/B* | K00370/K00371/K00374; K02567/K02568 | (as §5) | (as §5) | (as §5) | Molybdopterin oxidoreductase | NO₃⁻ → NO₂⁻ | Same nitrate-reducing enzymes as denitrification — the pathways diverge only at nitrite |
| Cytochrome c nitrite reductase (ammonia-forming) | *nrfA* | K03385 | PF02335 | IPR003321 | TIGR03152 | Pentaheme cytochrome c552 (CCNiR) | **NO₂⁻ → NH₄⁺** (committing step) | The defining DNRA marker — six-electron reduction of nitrite straight to ammonium. Its presence (vs *nirK/nirS*) flags N-retentive metabolism |
| NrfA accessory (membrane anchor) | *nrfH* ⚠️ | K15876 | PF13442 ⚠️ (NapC/NirT cyt c) | IPR011031 | — | Tetraheme cytochrome c (NapC/NirT) | NO₂⁻ → NH₄⁺ (electron transfer) | Quinol-oxidizing partner of NrfA |
| NADH-dependent nitrite reductase, large | *nirB* ⚠️ | K00362 | PF01077 + PF03460 ⚠️ + PF18267 | IPR017941 | TIGR02374 | NIR_SIR siroheme | NO₂⁻ → NH₄⁺ (cytoplasmic) | Cytoplasmic DNRA route. ⚠️ **PF01077 (NIR_SIR) is shared with assimilatory `nirA` and dissimilatory sulfite reductase `dsrAB`** — BLAST gating essential |
| NADH-dependent nitrite reductase, small | *nirD* ⚠️ | K00363 | PF01077 / PF03460 ⚠️ | IPR017941 | — | NIR_SIR ferredoxin half | NO₂⁻ → NH₄⁺ | Small subunit of NirBD |

---

## 7. Anammox — Anaerobic Ammonium Oxidation (NH₄⁺ + NO₂⁻ → N₂)

Planctomycete-specific (*Brocadia*, *Kuenenia*, *Jettenia*, *Scalindua*) pathway responsible for a large share of marine N loss. Couples ammonium and nitrite via the signature intermediate **hydrazine**. KEGG module **M00973**.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Hydrazine synthase α | *hzsA* ⚠️ | K20932 | ⚠️ multiheme cyt c (no dedicated Pfam) | — | — | Multiheme cytochrome c | NO + NH₄⁺ → N₂H₄ | **The definitive anammox phylomarker** (used as the field's diagnostic gene). Synthesizes hydrazine — chemistry unique to anammox. ⚠️ Detect via KO + BLAST (no stable Pfam; per-subunit custom HMMs were attempted but reverted — see config/manifests). ⚠️ **KEGG↔UniProt naming conflict:** K20932 ("hzsA") drives the call, but the ~809-aa α-subunit that UniProt labels HzsA is the one KEGG attaches to **K20934** — so the gene↔KO↔subunit-size mapping is scrambled vs the Kartal 2011 nomenclature. The subunits already score F1 1.00 via KO+BLAST and have no homology trap, so BLAST-only is the honest classification |
| Hydrazine synthase β | *hzsB* ⚠️ | K20933 | ⚠️ multiheme cyt c | — | — | Multiheme cytochrome c | N₂H₄ synthesis | Obligate subunit of the HZS heterotrimer |
| Hydrazine synthase γ | *hzsC* ⚠️ | K20934 | ⚠️ multiheme cyt c | — | — | Multiheme cytochrome c | N₂H₄ synthesis | Obligate subunit of the HZS heterotrimer |
| Hydrazine dehydrogenase / oxidoreductase | *hdh* (*hzo*) ⚠️ | K20935 | ⚠️ octaheme (HAO-like) | — | — | Octaheme cytochrome c | N₂H₄ → N₂ | Oxidizes hydrazine to dinitrogen, recovering electrons. ⚠️ HAO-family fold — BLAST-gate against anammox Hdh seeds |
| Nitrite reductase (NO supply) | *nirS* / *nirK* | K15864 / K00368 | (as §5) | (as §5) | (as §5) | (as §5) | NO₂⁻ → NO | Supplies the NO co-substrate for hydrazine synthesis; *Kuenenia* uses a *nirS*, others a *nirK*-like enzyme |

---

## 8. Assimilatory Nitrate / Nitrite Reduction (NO₃⁻ → NH₄⁺ for biosynthesis)

Reduction of nitrate/nitrite to ammonium for incorporation into biomass (not energy conservation). Cytoplasmic, regulated by N status. KEGG module **M00531**.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Ferredoxin-nitrate reductase | *narB* | K00367 | PF00384 + PF01568 | IPR006657 | — | Molybdopterin oxidoreductase | NO₃⁻ → NO₂⁻ (assimilatory) | Cyanobacterial/plant ferredoxin-dependent assimilatory nitrate reductase |
| Assimilatory nitrate reductase catalytic | *nasA* | K00372 | PF00384 + PF01568 | IPR006657 | — | Molybdopterin oxidoreductase | NO₃⁻ → NO₂⁻ (assimilatory) | Bacterial NADH assimilatory nitrate reductase (large/catalytic) |
| Assimilatory nitrate reductase, electron transfer | *nasB* | K00360 | PF07992 | IPR023753 | — | Pyr_redox (NADH oxidation) | NO₃⁻ → NO₂⁻ | Diaphorase partner of NasA |
| Nitrate reductase (NAD(P)H), eukaryotic | *NR / nia* | K10534 | PF00174 + PF00173 + PF00970 | IPR014756 | — | Mo-oxidoreductase + cyt b5 + FAD | NO₃⁻ → NO₂⁻ | Plant/fungal/algal assimilatory nitrate reductase |
| Ferredoxin-nitrite reductase | *nirA* ⚠️ | K00366 | PF01077 + PF03460 ⚠️ | IPR017941 | — | NIR_SIR siroheme | NO₂⁻ → NH₄⁺ (assimilatory) | Assimilatory nitrite reduction. ⚠️ Shares PF01077 with *nirB/nirD* (DNRA) and *dsrAB* — distinguish by ferredoxin vs NAD(P)H dependence + BLAST |
| Assimilatory nitrite reductase (NAD(P)H) | *nasD/nasE* | K17877, K00362/K00363-like, K26138/K26139, K00361 | PF01077 + PF03460 ⚠️ | IPR017941 | — | NIR_SIR siroheme | NO₂⁻ → NH₄⁺ | NAD(P)H-dependent assimilatory nitrite reductase |

---

## 9. Ammonia Assimilation (NH₄⁺ → glutamine / glutamate)

Incorporation of ammonium into the central metabolite glutamate — the gateway of N into all biosynthesis. GS/GOGAT (high-affinity) vs GDH (low-affinity) routes.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Glutamine synthetase | *glnA* | K01915 | PF00120 + PF03951 | IPR008147 | TIGR00653 | Gln-synt (GS) | NH₄⁺ + Glu → Gln | Primary high-affinity ammonium-assimilation enzyme; near-universal. The dominant route at low NH₄⁺ |
| Glutamate synthase (GOGAT) large | *gltB* | K00265 | PF01645 + PF04898 + PF01493 | IPR002489 | TIGR01317 | Glu_synthase (FMN/Fe-S) | Gln + 2-OG → 2 Glu | Completes the GS/GOGAT cycle; large catalytic subunit |
| Glutamate synthase (GOGAT) small | *gltD* | K00266 | PF07992 | IPR023753 | — | Pyr_redox (NAD(P)H) | Gln + 2-OG → 2 Glu | Electron-transfer subunit; obligate partner of GltB |
| Glutamate dehydrogenase | *gdhA* | K00261 (NAD(P)); K00262 (NADP) | PF00208 + PF02812 | IPR006095 | — | ELFV dehydrogenase | NH₄⁺ + 2-OG → Glu | Low-affinity assimilation at high NH₄⁺; reversible (also catabolic) |

---

## 10. Organic-N Mineralization / Ammonification (organic N → NH₄⁺)

Release of ammonium from organic nitrogen. **Urease** is the canonical, well-conserved marker (urea hydrolysis); broader amino-acid/amide deamination is taxonomically diffuse.

| Enzyme | Gene | KO | Pfam | InterPro | TIGRFAM/NCBIfam | Domain/Family | Step in cycle | Biogeochemical relevance |
|---|---|---|---|---|---|---|---|---|
| Urease α (catalytic) | *ureC* | K01428 | PF00449 + PF01979 | IPR006680 | TIGR01792 | Urease_alpha (Ni active site) | Urea → 2 NH₃ + CO₂ | The mineralization marker of choice; the catalytic Ni-binding subunit. Drives the largest single organic-N → NH₄⁺ flux in soils/sediments |
| Urease β | *ureB* | K01429 | PF00699 | IPR002019 | TIGR00192 | Urease_beta | Urea hydrolysis (structural) | Structural subunit of the (αβγ)₃ urease |
| Urease γ | *ureA* | K01430 | PF00547 | IPR002026 | TIGR00193 | Urease_gamma | Urea hydrolysis (structural) | Structural subunit |
| Urease accessory (Ni insertion) | *ureD/E/F/G* | K03187, K03188, K03189, K03190 | PF01774 / PF02814 / PF01730 / PF02492 | — | — | Urease maturation | Cofactor assembly | Required to activate apo-urease; their presence confirms a functional urease system |

---

## Obligatory Complexes — Must Find All Members Together

Mirrors **Info-ewaste.md** Part 1. These proteins are non-functional in isolation — multi-subunit catalytic assemblies where a single-subunit hit in a MAG is a red flag (or, for the ⚠️ traps, a likely mis-annotation).

| Complex | Components | Reaction | Why obligatory |
|---|---|---|---|
| **Nitrogenase (Mo)** | *nifH* + *nifD* + *nifK* | N₂ → NH₃ | NifH delivers electrons; NifD+NifK form the α₂β₂ FeMo-co catalytic core. None functions alone. *nifH*-only hits without *nifDK* (and ideally *nifENB*) indicate an incomplete/non-functional system |
| **Ammonia monooxygenase** | *amoA* + *amoB* + *amoC* | NH₃ → NH₂OH | AmoCAB is a single membrane Cu-monooxygenase complex. ⚠️ All three subunits are shared with methanotroph PmoCAB — the **set** must be co-validated against *amoA*-clade, not *pmoA*-clade, seeds |
| **Nitrite oxidoreductase** | *nxrA* + *nxrB* | NO₂⁻ → NO₃⁻ | α (catalytic, Mo) + β (Fe-S) form the minimal NXR. ⚠️ Both subunits share KOs with NarGH — co-occurrence plus NOB-clade BLAST is required to call NXR vs Nar |
| **Membrane nitrate reductase** | *narG* + *narH* + *narI* | NO₃⁻ → NO₂⁻ | Catalytic α + Fe-S β + membrane γ; γ anchors the complex and oxidizes quinol. Missing *narI* = no membrane electron input |
| **Periplasmic nitrate reductase** | *napA* + *napB* | NO₃⁻ → NO₂⁻ | Catalytic NapA depends on the di-heme NapB for electron delivery in the periplasm |
| **Nitric oxide reductase (cNOR)** | *norB* + *norC* | NO → N₂O | NorB (heme-Cu catalytic) + NorC (cytochrome c electron donor) form cNOR. (qNOR/*norZ* is a single-subunit alternative that does not require *norC*) |
| **Hydrazine synthase** | *hzsA* + *hzsB* + *hzsC* | NO + NH₄⁺ → N₂H₄ | The HZS heterotrimer; hydrazine synthesis requires all three subunits assembled. The obligate signature of anammox |
| **Glutamate synthase (GOGAT)** | *gltB* + *gltD* | Gln + 2-OG → 2 Glu | Large catalytic (GltB) + electron-transfer (GltD) subunits; non-functional apart |
| **Urease** | *ureA* + *ureB* + *ureC* (+ accessory *ureDEFG*) | urea → NH₃ | (αβγ)₃ assembly; accessory proteins insert the catalytic Ni. Structural subunits without UreC, or apo-enzyme without accessories, are inactive |

---

## Process-Completeness & Synergistic Interpretation

Mirrors **Info-ewaste.md** Part 2. Beyond single complexes, the ecological/biogeochemical meaning of a MAG comes from which *module* it completes. These are scored as completeness fractions across multiple genes.

| Pattern | Genes required | Interpretation / why it matters |
|---|---|---|
| **Complete denitrification** | (*narG/H/I* or *napA/B*) → (*nirK* or *nirS*) → (*norB/C*) → *nosZ* | A genome with all four steps fully respires NO₃⁻ to N₂. Most environmental denitrifiers are **partial** — the distribution of who has which step controls N-loss vs intermediate accumulation |
| **N₂O source vs sink** | has *nir*+*nor* but **lacks** *nosZ* (source) ↔ has *nosZ* (sink, esp. *nosZ*-II without upstream genes) | The single most climate-relevant call: organisms without *nosZ* leak the greenhouse gas N₂O; *nosZ*-II-only organisms are dedicated N₂O sinks |
| **Comammox completeness** | *amoCAB* + *hao* + *nxrAB* in **one** MAG (comammox *amoA*-clade) | Distinguishes a single-organism complete nitrifier (*Nitrospira*) from the classical two-step AOB/AOA + NOB division of labor |
| **Anammox completeness** | *hzsABC* + *hdh* + (*nirS* or *nirK*) | Confirms a functional anammox bacterium vs spurious *hzs* fragments; the marine/wastewater N-removal phenotype |
| **DNRA vs denitrification branch** | nitrate reductase + *nrfA* (DNRA, retains N) ↔ nitrate reductase + *nirK/S* (denitrification, loses N) | At the nitrite node, *nrfA* (→NH₄⁺) vs *nir* (→NO→gas) decides whether an ecosystem retains or loses fixed nitrogen — a major control on N availability |
| **Nitrifier-denitrifier coupling** | *amoCAB*/*hao* + *nirK* + *norB* (in AOB) | Some ammonia oxidizers also reduce their own nitrite to N₂O ("nitrifier denitrification") — a cryptic N₂O source flagged by *nir/nor* co-occurring with *amo* |
| **Full assimilatory branch** | (*narB* or *nasAB*) + (*nirA* or *nasDE*) + GS/GOGAT (*glnA*+*gltBD*) | Organism can build biomass N from nitrate end-to-end |

---

## Recommended Tiered Detection Strategy (for the MAG-mapping tool)

Following the **ewaste pipeline** architecture (Prodigal → HMMER → BLAST gating → complex/synergy scoring), with KO as the added primary anchor:

**Tier 1 — KO + well-behaved Pfam/TIGRFAM (direct calls).** Genes whose KO and HMM are clade-specific enough to call directly: *nifH* (PF00142/TIGR01287), *nifD/K*, *nirK* (TIGR02376), *nirS* (TIGR02375), *nrfA* (TIGR03152), *nosZ* (TIGR04244), *narI*, *napB*, *glnA* (TIGR00653), *gltB* (TIGR01317), *ureC* (TIGR01792).

**Tier 2 — KO + custom HMM + BLAST gate (the ⚠️ homology traps).** Require curated custom HMMs and BLAST confirmation against clade seeds:
- *nxrA/nxrB* vs *narG/narH* — identical KOs K00370/K00371; gate on NOB-clade vs denitrifier-clade seeds + operon context.
- *amoA* vs *pmoA* (and AOB vs AOA vs comammox *amoA*) — PF02461/PF12942; gate on clade.
- *nirB/nirD/nirA* vs *dsrAB* — shared PF01077 NIR_SIR fold; gate on NAD(P)H/ferredoxin specificity.
- *napA* vs other Mo-oxidoreductases — broad PF00384+PF01568; dsrA-style gate (per ewaste `napA`).

**Tier 3 — KO + BLAST-only (no stable Pfam).** Multiheme markers lacking clean HMMs: *hzsA/B/C*, *hdh*, *hao*. Detect by KO + BLAST against curated reference seeds (the lanmodulin/golB precedent in ewaste).

---

## Canonical Reference Organisms (for UniProt seed curation)

| Process | Model organism(s) | Notes |
|---|---|---|
| N₂ fixation | *Azotobacter vinelandii*, *Klebsiella pneumoniae* | Classic *nif* operon references |
| Ammonia ox. (AOB) | *Nitrosomonas europaea* | Canonical *amoCAB*/*hao* |
| Ammonia ox. (AOA) | *Nitrosopumilus maritimus*, *Nitrososphaera viennensis* | Archaeal *amoA* (PF12942) |
| Nitrite ox. (NOB) | *Nitrobacter winogradskyi*, *Nitrospira moscoviensis* | *nxrAB* seeds |
| Comammox | *Nitrospira inopinata* | Comammox *amoA* clade reference |
| Denitrification | *Paracoccus denitrificans*, *Pseudomonas stutzeri* | Full *nar/nap/nir/nor/nos* |
| DNRA | *Escherichia coli* (*nrfA*), *Wolinella succinogenes* | CCNiR references |
| Anammox | *Kuenenia stuttgartiensis*, *Brocadia* spp. | *hzsABC*, *hdh* |
| Assimilatory | *Bacillus subtilis* (*nasABDE*), *Synechocystis* (*narB/nirA*) | — |
| GS/GOGAT/GDH | *Escherichia coli* | *glnA*, *gltBD*, *gdhA* |
| Urease | *Klebsiella aerogenes*, *Helicobacter pylori* | Canonical urease structure |

---

## Sources

- **KEGG modules** (KO definitions verified against the REST API): [M00175 N fixation](https://www.genome.jp/module/M00175), [M00528 Nitrification](https://www.genome.jp/module/M00528), [M00529 Denitrification](https://www.genome.jp/module/M00529), [M00530 DNRA](https://www.genome.jp/module/M00530), [M00531 Assimilatory nitrate reduction](https://www.genome.jp/module/M00531), [M00804 Comammox](https://www.genome.jp/module/M00804), M00973 Anammox; pathway map [map00910 Nitrogen metabolism](https://www.genome.jp/pathway/map00910).
- **InterPro / Pfam**: [PF02461 AmoA/PmoA](https://www.ebi.ac.uk/interpro/entry/pfam/PF02461/), [PF12942 archaeal AmoA](https://www.ebi.ac.uk/interpro/entry/pfam/PF12942/), and per-gene families above.
- **NCBI NCBIfam / TIGRFAM** HMM accessions (TIGR01287 nifH, TIGR01580 narG/nxrA, TIGR04244 nosZ, TIGR03152 nrfA, etc.).
- **NCycDB** (Tu et al. 2019, *Bioinformatics*) — curated nitrogen-cycle gene-family database used to cross-check gene-family membership.
- Key reviews: Kuypers, Marchant & Kartal 2018 (*Nat Rev Microbiol*, "The microbial nitrogen-cycling network"); Daims et al. 2015 & van Kessel et al. 2015 (comammox); Kartal et al. 2011 (hydrazine synthase).

> **ID-verification status:** all KO numbers verified against the KEGG REST API (May 2026). Pfam/InterPro/TIGRFAM marked ⚠️ are either shared across reactions (homology traps) or lack a dedicated profile (multiheme cytochromes) and must be confirmed/curated with a custom HMM + BLAST gate during the tool-build phase — not relied on as a sole signature.
