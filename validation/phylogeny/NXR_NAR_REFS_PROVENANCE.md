# NxrA vs NarG type-labeled reference set — provenance & curation rationale

**Purpose.** An ORTHOGONAL, phylogeny-anchored ground truth to resolve the
identity/direction of a type-II DMSO-reductase-family Mo-bis-PGD α-subunit:
**NxrA** (nitrite oxidoreductase α, nitrite oxidizers / NOB, NO2⁻→NO3⁻) vs
**NarG** (membrane-bound respiratory nitrate reductase α, denitrifiers /
nitrate reducers, NO3⁻→NO2⁻). NxrA and NarG are paralogs that share KEGG KO
**K00370** and Pfam architecture — they are the same "trap" the ncycle pipeline
resolves with a NOB-clade custom HMM + a strict 60% BLAST gate + operon synteny.
A best-hit / placement call against this type-labeled set is therefore a
*non-circular* check on the pipeline's NxrA↔NarG verdict.

**Files**
- `nxr_nar_refs.faa` — 20 curated reference proteins, headers `>{acc}|{gene}|{type}|{organism}`.
- `nxr_nar_refs.tsv` — `accession  gene(nxrA|narG)  type(nxr|nar)  subtype  organism  source  notes`.
- `nxr_nar_refs.dmnd` — DIAMOND v2.2.0 database for the lightweight best-hit tier.
- Build scripts: `_query_nxrnar.py`, `_supp_nxr.py`, `_materialize_nxrnar.py`.

## Composition (verified-fetched 2026-06-13)

| gene | type | n |
|------|------|---|
| nxrA | nxr  | 12 |
| narG | nar  | 8  |
| **total** | | **20** |

Subtypes (NOB clades + canonical NarG):
Nitrospira-type 4 · Nitrobacter-type 3 · Nitrococcus-type 1 · Nitrospina-type 1 ·
Nitrolancea-type 1 · Nitrotoga-type 1 · anammox-type 1 · canonical_NarG 8.

All sequences fetched from UniProt REST, restricted to the expected α-subunit
length window (900–1350 aa) and protein/gene-name match. NarG entries are all
Swiss-Prot *reviewed*; most NxrA entries are *unreviewed* (TrEMBL) because NXR is
under-curated and frequently mis-annotated (see caveat 1). No accession invented.

### Clade coverage
- **NxrA (nitrite oxidizers):**
  - *Nitrobacter*-type (Alphaproteobacteria): *N. winogradskyi* (×2), *N. hamburgensis* — the
    cytoplasmic-facing NXR most similar to NarG (the hardest residual).
  - *Nitrococcus mobilis* (Gammaproteobacteria).
  - *Nitrospira*-type (periplasmic, divergent): *Ca. N. defluvii*, *N. japonica*,
    *Ca. N. nitrificans*, comammox *Ca. N. inopinata*.
  - *Nitrospina watsonii* (deep-branching marine NOB, periplasmic/divergent NXR).
  - *Nitrolancea hollandica* (Chloroflexi NOB).
  - *Ca. Nitrotoga* (Betaproteobacteria NOB).
  - *Ca. Brocadia fulgida* (anammox; the NXR/Nar-like enzyme that oxidizes NO2⁻→NO3⁻).
- **NarG (denitrifiers / nitrate reducers):** *E. coli*, *Pseudomonas aeruginosa*,
  *Bacillus subtilis*, *Thermus thermophilus*, *Mycobacterium tuberculosis*,
  *Halomonas maura*, *Salmonella* Typhimurium, *Geobacillus stearothermophilus* —
  Proteobacteria + Firmicutes + Actinobacteria + Deinococcus-Thermus spread.

## Superfamily relationships (the type-II DMSO reductase / Mo-bis-PGD family)

NxrA and NarG belong to the **type-II (membrane-bound) clade of the DMSO reductase
molybdoenzyme superfamily** (Mo-bis-molybdopterin-guanine-dinucleotide; also called
the complex-iron-sulfur-molybdoenzyme / CISM family). Their closest relatives in
this family — EbdA, DdhA, and especially the broad selenate/chlorate/polysulfide
reductases, plus the periplasmic NapA (type-I) — can also be pulled by a bare
K00370/Pfam search, which is why a *type-labeled* reference and a margin-based call
are needed rather than raw KO presence. NxrA and NarG cannot be separated below
~60% amino-acid identity (the basis of ncycle's 60% nxrA gate vs 45% narG gate).
The reference subtypes (Nitrobacter-type ≈ NarG-like; Nitrospira/Nitrospina-type =
periplasmic and more divergent) encode exactly the gradient of separability a
reviewer cares about.

## Circularity note (ncycle presence gate)

ncycle's nxrA presence BLAST gate (`config/targets.yaml`) uses
`Q71RT9, A0A0S4LQF4, A0A7S8FGD0` (nxrA) and `Q51075, A0A0S4LK29, M1KVL1` (nxrB).
On 2026-06-13, **Q71RT9 and A0A7S8FGD0 did not resolve via UniProt REST**
(demerged/obsolete); only **A0A0S4LQF4** (*Ca. Nitrospira nitrificans* nxrA1) is
live and **is included here, flagged `[ncycle-presence-ref]`** in the TSV. The nxrB
gate refs are a different subunit and are not part of this α-subunit set. The
direction/identity reference is otherwise fully independent of ncycle's presence
seeds (20 vs the 1 overlapping live nxrA seed), spanning every NOB clade + a broad
denitrifier NarG panel, so a NxrA↔NarG call from this set does not reduce to
ncycle's detection seeds.

## KNOWN CAVEATS (what a reviewer will raise)

1. **NXR is systematically mis-annotated as "nitrate reductase (quinone)".** Every
   *Nitrobacter* and *Nitrococcus* NXR α-subunit in UniProt is annotated
   `nitrate reductase (quinone), EC 1.7.5.1` — i.e. NarG — because of the paralogy.
   We label these `nxr` based on the organism being an obligate NOB with no
   denitrification pathway (Starkenburg et al. 2006/2008 genome analyses) and the
   gene's NXR operon context, NOT on the UniProt protein name. Each carries the
   `UniProt-annot='nitrate reductase(quinone)'` caveat in its TSV note. This IS the
   trap the reference exists to resolve — but it means the labels are curatorial
   (literature-derived), not taken from the database annotation.
2. **Nitrobacter NXR ≈ NarG (the hardest residual).** The cytoplasmic-facing
   *Nitrobacter*/*Nitrococcus* NXR is the most NarG-like (~55–60% id); scores
   overlap with denitrifier NarG. A best-hit call near the 60% boundary can flip.
   Report a margin/placement-support, and treat the Nitrobacter-type clade as the
   known low-confidence zone.
3. **Periplasmic / divergent NXR (Nitrospira, Nitrospina) loses BLAST identity.**
   The *Nitrospira*-type and *Nitrospina* NXR face the periplasm and are divergent
   (often <60% to canonical NarG and to Nitrobacter NXR). ncycle's own REPORT notes
   *Nitrospina* NXR sits at ~54% BLAST (below the 60% gate) and is recovered by the
   NOB-clade custom HMM — a sequence-only best-hit tier may under-call these. The
   set includes them so placement (not bare %id) can catch them.
4. **anammox NXR/Nar ambiguity.** The anammox enzyme (here *Brocadia*) catalyzes
   NO2⁻→NO3⁻ (oxidation, NXR-like physiology) but is phylogenetically Nar-like in
   some treatments. We label it `nxr` by physiology with an explicit note; a
   reviewer may prefer it as a separate "anammox-Nar" category.
5. **NarG paralogs / second nar operons.** Some genomes carry two NarG (narG/narZ)
   and other type-II Mo-enzymes (selenate, chlorate, polysulfide, ethylbenzene
   dehydrogenase). A query proteome can best-hit a non-Nar/non-Nxr type-II enzyme;
   gate on the NXR-vs-NarG margin AND exclude when neither clade hit is strong.
6. **NapA (periplasmic nitrate reductase) is a sibling, not in this set.** It is
   type-I and resolved separately in ncycle; do not call NapA hits as NarG. The
   length/name filter and DIAMOND target labels keep NapA out of the reference, but
   a raw search will still see it.
7. **Subtype labels are clade hypotheses, not flux measurements.** As with DsrAB,
   the call is "sequence-clade identity," not a measured direction.

## Obtaining a fuller phylogenetic-placement (EPA-ng / gappa) reference

The 20 isolates anchor each NOB clade + a denitrifier panel for a fast best-hit
tier. For a publication-grade placement tier:

- Build a curated type-II DMSO-reductase α-subunit alignment seeded from this set
  plus the canonical references in **Lücker et al. 2010** (*PNAS* 107:13479, the
  *Ca. Nitrospira defluvii* genome that defined the periplasmic Nitrospira-type
  NXR) and **Daims et al. 2015 / van Kessel et al. 2015** (comammox NXR), and the
  NarG references in the **MEROPS/InterPro IPR006468 (NarG) and IPR006657
  (Mo-bis-PGD)** families.
- **GraftM** ships nxrA/narG gpkgs; **PALADIN/PaPrBaG**-style or a custom
  `hmmalign → epa-ng → gappa examine assign` workflow against an NXR/NarG reference
  tree (rooted with NapA / formate dehydrogenase as outgroup) gives a placement
  with support values. Use the subtype taxonomy here as the gappa label set.
- For NOB taxonomy/markers see also **Pester et al. 2014** (Nitrospira nxrB
  phylomarker) — nxrB is more clade-resolved than nxrA, so a parallel nxrB
  reference can corroborate borderline Nitrobacter-type calls.

## Lightweight best-hit tier (this set)

```
diamond blastp --db nxr_nar_refs.dmnd --query <proteome.faa> \
  --outfmt 6 qseqid sseqid pident length bitscore stitle \
  --max-target-seqs 5 --evalue 1e-30 -k 5
```
Read the `type` token (`nxr`|`nar`) from the top hit's `sseqid`. Require the top
nxr-vs-nar bitscore margin to exceed a threshold and ideally top-N agreement;
treat Nitrobacter-type-vs-NarG ties (~55–60% id) as ambiguous and defer to the
placement tier + the NOB custom HMM (caveats 2–3).
