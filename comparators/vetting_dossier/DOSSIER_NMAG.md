# Nitrogen MAG realism study — candidate dossier (B7)

Mirrors scycle-pipeline's P3 MAG study (6 sulfur MAGs + 2 negatives, audited GO/NO-GO).
Purpose: exercise the nucleotide-mode **nxrA/narG operon-synteny resolver** (dormant on
isolate proteomes) on real N-cycle genomes spanning the trap. Candidates curated +
NCBI-verified by the selection specialist (2026-06-12). All accessions confirmed
downloadable nucleotide assemblies.

| # | Label | Accession | Organism | Phenotype | Environment / citation | Key N genes | Verdict |
|---|---|---|---|---|---|---|---|
| D1 | D1_Paracoccus_soil | GCF_000203895.1 | *Paracoccus denitrificans* PD1222 | denitrifier (complete) | soil model; PMID 17259312 | narGHI, napAB, nirS, norBC, nosZ | GO |
| D2 | D2_Thiobacillus_groundwater | GCF_000012745.1 | *Thiobacillus denitrificans* ATCC 25259 | denitrifier (S-driven) | groundwater; PMID 16452431 | narGHI, nirS, norBC, nosZ | GO — textbook narGHI operon |
| D3 | D3_Thioglobus_OMZ | GCF_001293165.1 | *Ca.* Thioglobus autotrophicus EF1 (SUP05) | denitrifier, N2O-emitter (no nosZ) | OMZ oxycline; PMID 26494660 | napAB, nirK, norBC (no nosZ) | GO — incomplete-denitrifier variant |
| N1 | N1_Nitrosopumilus_marine | GCF_000018465.1 | *Nitrosopumilus maritimus* SCM1 | AOA | marine; PMID 20421470 | amoABC, nirK | GO |
| N2 | N2_Nitrosomonas_AOB | GCF_000009145.1 | *Nitrosomonas europaea* ATCC 19718 | AOB | soil/wastewater; PMID 14583071 | amoABC, hao, nirK, norB | GO |
| N3 | N3_Nitrospina_NOB | GCF_000341545.2 | *Nitrospina gracilis* 3/211 | NOB | marine; PMID 23439773 | nxrAB | GO-WITH-CAVEAT — scaffold N50 57 kb, nxr split risk |
| N4 | N4_Nitrospira_comammox | GCF_001458695.1 | *Ca.* Nitrospira inopinata ENR4 | comammox | oil-well biofilm; PMID 26610024 | amoABC, hao, nxrAB | GO — prime amo+nxr dual-operon synteny test |
| A1 | A1_Brocadia_reactor | GCA_017347445.1 | *Ca.* Brocadia pituitae | anammox | anammox bioreactor MAG | hzsABC, hdh, nxr-type | GO — genuine reactor MAG |
| NC1 | Neg1_Bacteroides_gut | GCF_014131755.1 | *Bacteroides thetaiotaomicron* | negative (clean) | human gut isolate | none dissimilatory | GO — clean negative |
| NC2 | Neg2_Ecoli_K12 | GCF_000005845.2 | *Escherichia coli* K-12 MG1655 | negative (assimilatory/anaerobic only) | lab isolate | narG/narZ present, no nxr/nir-denit/nosZ | GO — "narG-present-not-denitrifier" stress negative |

**Phenotype coverage:** 3 denitrifiers (complete / S-driven / N2O-emitter) + AOA + AOB + NOB +
comammox + anammox + 2 negatives, across soil/groundwater/marine-OMZ/oil-well/wastewater.

**Best nxrA/narG synteny exercisers:** N4 comammox (amoABC **and** nxrAB — hardest), N3
Nitrospina + A1 Brocadia (nxr+, no narGHI), vs D1 Paracoccus + D2 Thiobacillus (canonical
narGHI with narI proximity) — true nxr-vs-narG positive pairs on both sides of the K00370 trap.
NC2 *E. coli* is the narG-present-but-not-a-denitrifier negative stress case.

## Notes / corrections from the specialist audit
- GCF_001458695.1 = *N. inopinata* ENR4 (NOT japonica; japonica = GCF_900169565.1, available as
  an alternate second comammox if wanted).
- Brocadia primary-publication PMID not pinned from the assembly record; alternate verified
  anammox reactor MAGs: *Ca.* Kuenenia stuttgartiensis GCA_011066545.1 / GCA_900232105.1.
- Several functional genomes are closed enrichment/reference genomes (canonical N-cycle
  phenotypes exist mainly as closed genomes); nucleotide input still exercises the synteny
  resolver. Adding 1–2 genuinely fragmented denitrifier MAGs would strengthen the
  fragmentation-realism angle (refinement; not blocking).
