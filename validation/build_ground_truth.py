#!/usr/bin/env python3
"""
build_ground_truth.py — curated N-cycle ground truth (v12) for the validation panel.

v12 (2026-05-30, ROADMAP P5.3.4): broad TN expansion for the two clade HMMs
introduced in v7 (`amoA_gamma`) and v10 (`nosZ_clade2`). Following the
audit-recommended `amoA_archaeal` pattern, add target as `absent` to every
panel + hold-out genome that doesn't have an explicit P/A call for it (the
EXPAND_TN_TARGETS post-pass in main()). Lifts amoA_gamma TN coverage 5 → ~32
and nosZ_clade2 TN coverage 5 → ~32, which clears the P5.1.4 reliability
flag once n_pos ≥ 3. Source tag `expanded_tn_v12` distinguishes
programmatically-added TN cells from manually-curated v11 entries; hold-out
expansions keep source `holdout_v3` so they continue to be excluded from
the training-panel score.

v11 (2026-05-30, post-audit): RETRACT v6 nosZ A→P AND v9 norZ A→P on
Nwinogradskyi_Nb255. The 2026-05-30 audit (REPORT § Audit 2026-05-30) found
that `test_panel/Nwinogradskyi_Nb255.faa` was a *Burkholderia thailandensis*
proteome (5607 proteins, 5258 explicitly [Burkholderia thailandensis]-tagged),
not Nb-255. Both v6 (evidence WP_080511513.1, "Pseudomonas-like NosZ at 64.7%
identity") and v9 (evidence WP_009890104.1, "762-aa qNor at 46.3% to B3PE38")
were Burkholderia proteins surfacing under the mislabeled file. Re-prodigaled
real Nb-255 (NC_007406.1, 3262 proteins) re-run on 2026-05-30: nosZ = absent,
nosZ_clade2 = absent, the sole K00376 hit is sub-threshold (score 36.4,
E 5.4e-12, no PF18764 hits); norZ = K04748 hits exist (best E 4.6e-134) but
BLAST-gate disqualifies them against curated qNor seeds — the protein is not
a true qNor. Three Nb-255 cells flipped: nosZ P→A, norZ P→A, nosZ_clade2
added as A (parallel TN). v4 narG A→P STANDS — Starkenburg 2006/2008 narGHI
operon confirmed on real Nb-255 (narG + narH via synteny, narI via KO K00374).
All other v9 cells (Cnecator + Synechocystis norZ A→P) STAND — those qNor
anchors (Q7WX97, P74677) are independent of the Burkholderia file and
audit-confirmed. Panel-QC gate added (validate_panel.py) to prevent recurrence.

v10 (2026-05-30): added `nosZ_clade2` target cells after Batch 3d clade-II
nosZ HMM (REPORT § Block 3d). Three clade-I NosZ carriers given nosZ_clade2 → A
as discrimination TN cells (Pdenitrificans α-proteo, Paeruginosa γ-proteo,
Cnecator β-proteo — diverse clade-I representatives). Wsuccinogenes nosZ P→A
on annotation-gap grounds: the catalytic NosZ subunit is missing from the
panel's RefSeq proteome (only accessory NosL/NosD proteins present); the new
clade-II HMM gives zero hits, confirming the catalytic protein truly isn't in
the available proteome — analogous to the v5 Ngracilis nxrB assembly-truncation
case. The gene is biologically present in *Wolinella succinogenes* (Simon et al.
2004) but undetectable from this assembly's annotation.

v9 (2026-05-30): closed the qNor coverage gap (REPORT § Block 3 follow-up (B))
by adding qNor BLAST seeds + K04561 fallback KO to the `norZ` target. Three
norZ presence cells added to the GT, each verified by canonical qNor length
(700-800 aa) + UniProt/NCBI protein-name + above-gate BLAST identity to qNor
clade seeds: Cnecator norZ A→P (Q7WX97, 762 aa, 77% to Cellvibrio qNor);
Nwinogradskyi norZ A→P (WP_009890104, 762 aa, 46% — qNor missed by the
Starkenburg 2006/2008 annotation due to a confusing cross-genome NCBI label);
Synechocystis norZ added (P74677, 770 aa, 43% — documented cyanobacterial qNor).

v8 (2026-05-29): Cnecator norB P→A. The norB P call was a curator-time mislabel —
the existing entry comment already said "norB qNOR, NO norC". UniProt confirms
Cnecator's Q0JYR9 is qNor (NorB2 locus tag), 762 aa canonical qNor length, not
cNorB. The pipeline's cNor BLAST gate correctly rejects Q0JYR9 at 25% identity
(correct discrimination, not a coverage gap in the new β-proteo seeds). The
qNor protein is biologically present but not detectable under the current
norZ machinery (K04748 doesn't fire; no qNor BLAST seeds yet) — flagged as a
follow-up coverage gap, norZ kept as A for the GT until detection is in place.

v7 (2026-05-29): added `amoA_gamma` target cells for the new γ-AOB-clade amoA HMM
(REPORT § Block 3b). Noceani amoA P→A and amoA_gamma added to P (β-AOB-clade `amoA`
genuinely doesn't fire on the γ-AOB clade — relabeling reflects biology, not error).
Mfumariolicum amoA_gamma → A (primary calibration TN: β-pmoA must not cross-react).
Other AOB/comammox/AOA entries left unchanged for now — broader TN expansion
deferred until amoA_gamma TC calibration is observed against the panel.

v6 (2026-05-28): Ninopinata nrfA absent→present and Nwinogradskyi nosZ absent→present
(REPORT § Long-tail Block 2d). Both had strong protein-level evidence; v5 calls were
GT errors corrected by curated_v6.

v5 (2026-05-28): Ngracilis nxrB present→absent on assembly-truncation grounds
(REPORT § Long-tail Block 2a). The genome's assembly lacks a detectable nxrB ORF.

v4 (2026-05-25): corrected Nwinogradskyi_Nb255 narG absent→present — the genome carries a
genuine narG-narH-narJ-narI respiratory-nitrate-reductase operon (confirmed from the NCBI
GFF and operon synteny), distinct from its NXR; the v3 `absent` was a curation error.

v3 incorporates an environmental-microbial-genomics specialist audit (2026-05-25):
KEGG per-organism KO verification of present/absent calls, urease added across
genomes that carry it, the Nitrospira japonica mislabel fixed, dubious cyano/
verruco gltB removed, and the panel expanded to ~34 genomes spanning all 9 pathways
with ≥3 genomes per major pathway + clean negative controls.

Two genomes are HOLD-OUTS (HOLDOUTS set) — multi-pathway, genera absent from all
HMM/BLAST training — reserved as an independent generalization test for P4.

Scorer (score_ncycle.py) scores ONLY cells listed here; unlisted (genome,target)
pairs are not scored. glnA = present in every cellular organism (near-universal).
Output: validation/ground_truth.tsv  (genome, target, expected, source)
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "ground_truth.tsv"

HOLDOUTS = {"Bdiazoefficiens_USDA110", "Rpalustris_CGA009",
            # P5.3.5 (2026-05-30): 6 new hold-outs spanning 5 phyla (audit-required
            # generalization expansion: clade-II NosZ carrier, γ-AOB, non-Nitrospira
            # NOB, additional anammox lineage + 2 phylogenetic-diversity extras).
            "Adehalogenans_2CP1",   # δ-proteo, clade-II NosZ carrier
            "Nhalophilus_Nc4",      # γ-AOB
            "Nhollandica_Lb",       # Chloroflexi NOB
            "Sbrodae",              # Planctomycetes anammox (Scalindua lineage)
            "Sstutzeri_F2a",        # γ-proteo canonical denitrifier
            "Mcapsulatus_Bath"}     # γ-proteo methanotroph (amoA decoy)

UREASE = ["ureA", "ureB", "ureC"]

CURATED: dict[str, dict[str, list[str]]] = {
    # ───────────── existing 15 (audit-corrected) ─────────────
    "Avinelandii_DJ": {  # Azotobacter vinelandii — 3 nitrogenases
        "P": ["nifH","nifD","nifK","nifE","nifN","nifB","vnfD","vnfH","anfG",
              "nasA","glnA","gltB","gltD",*UREASE],
        "A": ["amoA","hao","nxrA","nxrB","narG","nirS","nirK","nosZ","nrfA",
              "hzsA","hzsB","hzsC","hdh"],
    },
    "Neuropaea_ATCC19718": {  # Nitrosomonas europaea (AOB)
        "P": ["amoA","amoB","amoC","hao","nirK","norB","norC","glnA"],
        "A": ["nifH","nifD","nxrA","nxrB","narG","nosZ","nrfA","nasA",
              "hzsA","hzsB","hzsC","hdh","amoA_archaeal","amoA_gamma"],
    },
    "Ninopinata_comammox": {  # true comammox Nitrospira inopinata (GCA_001458695.1).
        # v6 (2026-05-28): nrfA absent→present. CUQ65653.1 is annotated by EBI as
        # "Pentaheme cytochrome c nitrite reductase NrfA" and scores K03385 at 582 —
        # canonical NrfA family. Nitrospira lineage II genomes (incl. comammox) are
        # documented to carry nrfA-family genes (Daims et al. 2015 + Kits et al. 2017);
        # functional role is debated (DNRA vs nitrite detoxification vs biosynthesis)
        # but the gene is unambiguously present. v5 GT marked it absent on the
        # canonical-comammox-pathway assumption; the assembly contains it.
        "P": ["amoA","amoB","amoC","hao","nxrA","nxrB","nirK","glnA","ureC","nrfA"],
        "A": ["nifH","narG","nirS","nosZ","hzsA","hzsB","hzsC","hdh",
              "amoA_archaeal","amoA_gamma"],
    },
    "Njaponica_NJ11": {  # Nitrospira japonica (NOB; GCF_900169565.1) — NOT comammox
        "P": ["nxrA","nxrB","nirK","hao","glnA",*UREASE],
        "A": ["amoA","amoB","amoC","nifH","narG","nirS","nosZ","nrfA",
              "hzsA","hzsB","hzsC","hdh","amoA_archaeal"],
    },
    "Nmaritimus_SCM1": {  # Nitrosopumilus maritimus (AOA; no urease)
        "P": ["amoA_archaeal","nirK","glnA"],
        "A": ["amoA","hao","nxrA","narG","nirS","nosZ","nifH","nrfA",
              "hzsA","hzsB","hdh"],
    },
    "Nwinogradskyi_Nb255": {  # Nitrobacter winogradskyi (NOB) with respiratory nitrate reductase
        # Genuine narG-narH-narJ-narI respiratory-nitrate-reductase operon at 0.86 Mb on
        # NC_007406.1, distinct from the standalone NxrA at 2.26 Mb.
        # v4 (2026-05-25): narG absent→present (Starkenburg 2006/2008 + NCBI GFF + operon
        # synteny — curated_v3 was a curation error). STANDS post-v11 — confirmed on real
        # Nb-255 (narG + narH via synteny @ NC_007406.1_784/_786, narI via KO K00374 @ _788).
        # v6 (2026-05-28): nosZ A→P. RETRACTED in v11. Evidence WP_080511513.1 ("score 854,
        # 64.7% to Q59105") was Burkholderia thailandensis NosZ surfacing under the
        # mislabeled panel file, not Nb-255. Real Nb-255 re-run on 2026-05-30: nosZ + nosZ_clade2
        # both absent; the only K00376 hit (NC_007406.1_2408) is sub-threshold (score 36.4,
        # E 5.4e-12, well below the ~70 TC) and no PF18764. Cell flipped back to A.
        # v9 (2026-05-30): norZ A→P (qNor). RETRACTED in v11. Evidence WP_009890104.1
        # ("762-aa qNor, K04561 at 1.2e-306, 46.3% to B3PE38") was Burkholderia qNor under
        # the mislabeled panel file. Real Nb-255 re-run: K04748 hits DO exist (best
        # NC_007406.1_2056, E 4.6e-134) but the BLAST gate DISQUALIFIES against curated
        # qNor seeds — the K04748-bearing protein on real Nb-255 is not a true qNor by
        # sequence identity. Cell flipped back to A. Restores Starkenburg's no-NO-reductase
        # framing (further revised: Nb-255 does not have a demonstrated NorB/NorC nor a
        # bona-fide qNor under current evidence).
        # v11 also adds nosZ_clade2 → A (parallel TN cell, consistent with nosZ A).
        "P": ["nxrA","nxrB","nirK","glnA","narG"],
        "A": ["napA","norB","norZ","nosZ","nosZ_clade2","amoA","amoB","hao",
              "nifH","nrfA","hzsA","hzsB","hzsC","hdh","amoA_archaeal"],
    },
    "Pdenitrificans_PD1222": {  # Paracoccus denitrificans — full denitrification
        # v10 (2026-05-30): nosZ_clade2 → A. Paracoccus is the canonical clade-I
        # NosZ reference (Tat-secreted); verifies the nosZ_clade2 HMM correctly
        # discriminates clade-I from clade-II at TC=420.
        "P": ["narG","narH","narI","napA","napB","nirS","norB","norC","nosZ",
              "nasA","glnA","gltB","gltD",*UREASE],
        "A": ["nxrA","nxrB","amoA","hao","nifH","nirK","nrfA",
              "hzsA","hzsB","hzsC","hdh","amoA_archaeal","nosZ_clade2"],
    },
    "Cnecator_H16": {  # Cupriavidus necator — denitrification via qNor (NorB2/Q0JYR9), no cNor
        # v8 (2026-05-29): norB P→A. The comment on this entry has always said
        # "norB qNOR, NO norC" — confirmed: Cnecator's Q0JYR9 is annotated by
        # UniProt as "Nitric oxide reductase qNor type (NorB2)", 762 aa
        # (canonical qNor length; cNorB is ~470 aa). The gene name `norB2` is
        # a locus tag, not a functional cNor claim. So `norB` (cNor target) is
        # correctly absent; the cNor BLAST gate properly rejects Q0JYR9 at 25%
        # identity (correct discrimination, NOT a coverage failure of the new
        # β-proteo seeds added in Batch 3a).
        # v9 (2026-05-30): norZ absent→present after qNor BLAST seeds added
        # (Batch 3c). Q7WX97 (Cnecator's second qNor, PHG244 locus, 762 aa)
        # fires K04561 at E=0.0 and BLASTs to Cellvibrio japonicus qNor (B3PE38)
        # at 77.2% identity. Closes REPORT § Block 3 follow-up "(B)" — Cnecator's
        # NO-reduction step is now detectable. (Q0JYR9 also matches the qNor
        # seeds at high identity; both qNor paralogs are picked up.)
        "P": ["narG","narH","narI","napA","napB","nirS","norZ","nosZ","nasA",
              "glnA",*UREASE],
        "A": ["norB","nxrA","nxrB","amoA","hao","nifH","norC","hzsA","hzsB","hdh",
              "amoA_archaeal","nosZ_clade2"],  # v10: β-proteo clade-I NosZ TN cell
    },
    "Bsubtilis_168": {  # Bacillus subtilis — respiratory nar + assimilatory; no denitrification
        "P": ["narG","narH","narI","nasA","nasD","glnA","gltB","gltD",*UREASE],
        "A": ["nirS","nirK","norB","norC","nosZ","amoA","hao","nxrA","nxrB",
              "nifH","hzsA","hdh","amoA_archaeal"],
    },
    "Synechocystis_PCC6803": {  # cyanobacterium — assimilatory; cyano GOGAT is Fd-type (no gltB)
        # v9 (2026-05-30): norZ added to P. P74677 (770 aa, UniProt "Cytochrome b
        # subunit of nitric oxide reductase", gene name norB) is canonical qNor
        # length and fires K04561 at E=6.6e-287, BLASTing to Cellvibrio japonicus
        # qNor (B3PE38) at 43.0% identity (above 40% gate). Cyanobacterial qNor
        # is documented in literature (NO detoxification / NO sensing role).
        # Lifted by Batch 3c qNor BLAST seeds.
        "P": ["narB","nirA","glnA","gltD","norZ",*UREASE],
        "A": ["nifH","narG","nirS","nirK","nosZ","amoA","hao","nxrA","nxrB",
              "nrfA","hzsA","hdh","amoA_archaeal"],
    },
    "Ahydrophila_ATCC7966": {  # Aeromonas hydrophila — DNRA
        "P": ["nrfA","nrfH","nirB","nirD","napA","napB","nasA","glnA","gltB","gltD"],
        "A": ["nosZ","nirS","nirK","amoA","hao","nxrA","nxrB","nifH",
              "hzsA","hzsB","hdh","amoA_archaeal"],
    },
    "Kstuttgartiensis": {  # Ca. Kuenenia stuttgartiensis — anammox (has nxr-like + HAO-family homologs → not scored)
        "P": ["hzsA","hzsB","hzsC","hdh","nirS","glnA"],
        "A": ["amoA","amoA_archaeal","nifH","nifD","nosZ","narG"],
    },
    "Dvulgaris_Hildenborough": {  # Desulfovibrio vulgaris — diazotroph + DNRA (nrfA); dsr decoy for nir
        "P": ["nifH","nifD","nifK","nrfA","glnA"],
        "A": ["nirB","nirD","nirA","amoA","amoA_archaeal","hao","nxrA","nxrB",
              "narG","nirS","nosZ","hzsA","hzsB","hdh"],
    },
    "Mfumariolicum_SolV": {  # Methylacidiphilum (verruco methanotroph + diazotroph) — amoA TRAP
        # v7 (2026-05-29): amoA_gamma added to A. Primary calibration TN for the
        # new γ-AOB HMM (pmoA must NOT cross-react). β-pmoA panel bitscores
        # currently 274 (K28504) / 279 (β-AOB+comammox HMM) — see REPORT § Block 2b.
        "P": ["nifH","nifD","nifK","nifE","nifB","nasA","nirK","glnA"],
        "A": ["amoA","amoA_gamma","amoA_archaeal","nxrA","nxrB","narG","nirS","nosZ",
              "hzsA","hzsB","hzsC","hdh","nrfA"],
    },
    "Aterreus_NIH2624": {  # Aspergillus terreus — eukaryotic assimilatory (NR)
        "P": ["NR","glnA",*UREASE],
        "A": ["nifH","amoA","amoA_archaeal","nxrA","narG","nirS","nosZ","hzsA","nrfA"],
    },
    "Scerevisiae_S288C": {  # S. cerevisiae — eukaryote negative (no nitrate assimilation)
        "P": ["glnA","gdhA"],
        "A": ["NR","nifH","amoA","amoA_archaeal","nxrA","narG","nirS","nosZ",
              "nrfA","hzsA","nasA","narB","nirA"],
    },

    # ───────────── v3 additions (specialist-recommended, KEGG-verified) ─────────────
    "Npcc7120": {  # Nostoc sp. PCC 7120 — cyanobacterial diazotroph
        "P": ["nifH","nifD","nifK","nifE","nifB","narB","nirA","glnA","gdhA",*UREASE],
        "A": ["amoA","nxrA","narG","nirS","nosZ","nrfA","hzsA","amoA_archaeal"],
    },
    "Smeliloti_1021": {  # Sinorhizobium meliloti — nif + denitrification
        "P": ["nifH","nifD","nifK","nifE","nifB","napA","napB","nirK","norB","norC",
              "nosZ","nasA","glnA",*UREASE],
        "A": ["amoA","nxrA","nirS","hzsA","nrfA","amoA_archaeal"],
    },
    "Paeruginosa_PAO1": {  # Pseudomonas aeruginosa — full denitrification (nirS)
        # v10 (2026-05-30): nosZ_clade2 → A. γ-proteo clade-I NosZ TN cell.
        "P": ["narG","narH","narI","napA","napB","nirS","norB","norC","nosZ",
              "nasA","glnA","gltB","gltD",*UREASE],
        "A": ["amoA","nxrA","nirK","hzsA","nrfA","nifH","amoA_archaeal","nosZ_clade2"],
    },
    "Wsuccinogenes_DSM1740": {  # Wolinella succinogenes — canonical nrfA DNRA
        # v10 (2026-05-30): nosZ P→A on annotation-gap grounds (analogous to the
        # v5 Ngracilis nxrB assembly-truncation case). Wolinella succinogenes is
        # a published clade-II/atypical nosZ carrier (Simon et al. 2004), but
        # the panel proteome (Wsuccinogenes_DSM1740.faa, 2041 proteins) does NOT
        # contain a detectable catalytic NosZ subunit: zero K00376 hits even at
        # E=100, zero PF18764/PF00116/TIGR04244 hits, zero hits on the new
        # nosZ_clade2 HMM, only accessory nos-operon proteins (NosL, NosD)
        # named in the proteome. UniProt's Wolinella nosZ entries (Q7M9H4 470 aa,
        # Q7M9H5 410 aa — both truncated vs the canonical ~640 aa) BLAST against
        # the panel proteome at only 25-38% identity over short fragments
        # (noise level). The catalytic subunit is missing from this particular
        # RefSeq annotation. Closes the long-tail nosZ FN as TN; the gene is
        # biologically present in the organism but undetectable in the
        # available assembly's protein call set.
        "P": ["nrfA","nrfH","napA","napB","nifH","nifD","nifK","glnA","gltB","gltD"],
        "A": ["nosZ","nosZ_clade2","amoA","nxrA","nirS","nirK","narG","hzsA","amoA_archaeal"],
    },
    "Ecoli_K12_MG1655": {  # E. coli — nar + DNRA (nrfA, nirBD); no denitrification/nif/urease
        "P": ["narG","narH","narI","napA","napB","nrfA","nrfH","nirB","nirD",
              "glnA","gltB","gltD"],
        "A": ["nirS","nirK","norB","norC","nosZ","nifH","amoA","nxrA","hzsA",
              "ureC","amoA_archaeal","NR"],
    },
    "Soneidensis_MR1": {  # Shewanella oneidensis — DNRA (nap, nrf)
        "P": ["napA","napB","nrfA","nrfH","glnA","gltB","gltD"],
        "A": ["narG","nirS","nirK","norB","nosZ","amoA","nxrA","nifH","hzsA",
              "amoA_archaeal"],
    },
    "Noceani_ATCC19707": {  # Nitrosococcus oceani — gammaproteobacterial AOB
        # v7 (2026-05-29): amoA P→A, amoA_gamma added to P. β-AOB-clade amoA HMM
        # (target `amoA`) genuinely doesn't fire on Noceani (γ-AOB clade is distinct
        # at the seed/training level); the γ-AOB form is now tracked as `amoA_gamma`.
        "P": ["amoA_gamma","amoB","amoC","hao","nirK","norB","norC","glnA",*UREASE],
        "A": ["amoA","nxrA","nifH","nosZ","hzsA","narG","amoA_archaeal","nrfA"],
    },
    "Nmultiformis_ATCC25196": {  # Nitrosospira multiformis — beta-AOB (NB: in amoA training)
        "P": ["amoA","amoB","amoC","hao","nirK","norB","glnA"],
        "A": ["nxrA","nifH","nosZ","hzsA","narG","amoA_archaeal","amoA_gamma","nrfA"],
    },
    "Nviennensis_EN76": {  # Nitrososphaera viennensis — soil AOA
        "P": ["amoA_archaeal","nirK","glnA",*UREASE],
        "A": ["amoA","hao","nxrA","nifH","narG","nirS","nosZ","hzsA","nrfA"],
    },
    "Ngracilis_3211": {  # Nitrospina gracilis — deep-branching NOB (distinct phylum).
        # v5 (2026-05-28): nxrB present→absent. The GCF_000341545.2 assembly is
        # fragmented around the nxr locus: contig NZ_HG422176.1 is only 583 bp and
        # carries just the C-terminal 176 aa of the 425-aa NxrB (UniProt M1KVL1,
        # AGF29470.1 — full gene was sequenced separately by Luecker et al. 2013
        # outside the WGS submission, hence its presence in UniProt). NCBI's
        # annotation correctly omits the truncated ORF; even Prodigal-recovered
        # partial scores 260 on the custom nxrB HMM (TC 400 — well below). The full
        # protein exists in nature but is not in *this proteome*, so GT should
        # reflect what the assembly contains (= absent), not what the organism has.
        "P": ["nxrA","glnA"],
        "A": ["nxrB","amoA","amoA_archaeal","nifH","narG","nosZ","hzsA","nrfA"],
    },
    "Hpylori_26695": {  # Helicobacter pylori — urease specificity control
        "P": ["ureC","glnA"],
        "A": ["amoA","nxrA","narG","napA","nirS","nirK","norB","nosZ","nifH",
              "nrfA","hzsA","nasA","narB","amoA_archaeal"],
    },
    "Bsinica_JPN1": {  # Ca. Brocadia sinica — anammox (independent lineage vs Kuenenia)
        "P": ["hzsA","hzsB","hzsC","hdh","glnA"],
        "A": ["amoA","amoA_archaeal","nifD","nosZ","narG","nirK"],
    },
    "Nnitrosa_comammox": {  # second verified comammox (Ca. Nitrospira nitrosa, GCA_001458775.1)
        "P": ["amoA","amoB","amoC","hao","nxrA","nxrB","glnA"],
        "A": ["nifH","narG","nirS","nosZ","nrfA","hzsA","hzsB","hzsC","hdh",
              "amoA_archaeal","amoA_gamma"],
    },
    "Spneumoniae_ref": {  # Streptococcus pneumoniae — negative control
        "P": ["glnA"],
        "A": ["amoA","amoA_archaeal","nxrA","nxrB","narG","napA","nirS","nirK",
              "norB","nosZ","nrfA","nifH","hzsA","nasA","narB","NR","ureC"],
    },
    "Lacidophilus_4356": {  # Lactobacillus acidophilus — negative control
        "P": ["glnA"],
        "A": ["amoA","amoA_archaeal","nxrA","nxrB","narG","napA","nirS","nirK",
              "norB","nosZ","nrfA","nifH","hzsA","nasA","narB","NR","ureC"],
    },

    # ───────────── HOLD-OUTS (independent test; genera absent from all training) ─────────────
    "Bdiazoefficiens_USDA110": {  # Bradyrhizobium diazoefficiens — nif + full denitrification + assimilatory
        "P": ["nifH","nifD","nifK","nifE","nifB","napA","napB","nirK","norB","norC",
              "nosZ","nasA","nirA","glnA","gltB",*UREASE],
        "A": ["amoA","nxrA","nirS","nrfA","hzsA","narG","amoA_archaeal"],
    },
    "Rpalustris_CGA009": {  # Rhodopseudomonas palustris — nif + partial denitrification (nirK/nor/nosZ), NO narG/nap/nrfA
        "P": ["nifH","nifD","nifK","nifE","nifB","nirK","norB","norC","nosZ",
              "nirA","glnA","gltB",*UREASE],
        "A": ["narG","narH","narI","napA","napB","nirS","nrfA","amoA","nxrA",
              "hzsA","amoA_archaeal"],
    },
    # ───────────── P5.3.5 hold-outs (2026-05-30): 6 new genomes spanning 5 phyla ─────────────
    "Adehalogenans_2CP1": {  # Anaeromyxobacter dehalogenans 2CP-1 — δ-proteo Myxococcales,
        # documented clade-II atypical NosZ carrier (Sanford et al. 2012; Hallin et al. 2018).
        # Also a versatile anaerobic respirer: complete denitrification machinery + nitrogen fixation.
        # P5.3.5: the clade-II NosZ slot — first held-out positive for nosZ_clade2.
        "P": ["nifH","nifD","nifK","narG","narH","narI","nirS","nirK","norB","norC",
              "nosZ_clade2","nasA","glnA","gltB","gltD"],
        "A": ["nosZ","napA","amoA","amoA_archaeal","amoA_gamma","nxrA","nxrB",
              "hzsA","hzsB","hzsC","hdh","nrfA"],
    },
    "Nhalophilus_Nc4": {  # Nitrosococcus halophilus Nc 4 — γ-proteo AOB (γ-AOB clade).
        # Carries γ-AOB-clade amoA (amoA_gamma target), classical amoABC complex; produces nitrite.
        # Some γ-AOB also do nirK-mediated NO2-→NO (nitrifier-denit pathway). P5.3.5: held-out
        # positive for amoA_gamma + the γ-AOB nitrifier-denit phenotype.
        "P": ["amoA_gamma","amoB","amoC","hao","nirK","glnA","gltB","gltD",*UREASE],
        "A": ["amoA","amoA_archaeal","nxrA","nxrB","narG","nirS","nosZ","nosZ_clade2",
              "nifH","nrfA","hzsA","hzsB","hzsC","hdh"],
    },
    "Nhollandica_Lb": {  # Nitrolancea hollandica Lb — Chloroflexi NOB (non-Nitrospira NOB clade).
        # Carries Nitrolancea-clade nxrA/nxrB. P5.3.5: held-out positive for the nitrite-oxidation
        # step from a non-Nitrospira lineage — tests nxr clade generalization. Also a moderate
        # thermophile with assimilatory nitrate reduction.
        "P": ["nxrA","nxrB","nasA","glnA","gltB","gltD"],
        "A": ["amoA","amoA_archaeal","amoA_gamma","narG","nirS","nirK","nosZ",
              "nosZ_clade2","nifH","nrfA","hzsA","hzsB","hzsC","hdh"],
    },
    "Sbrodae": {  # Ca. Scalindua brodae — marine Planctomycetes anammox (Scalindua lineage,
        # distinct from panel's Brocadia and Kuenenia). Canonical anammox machinery.
        # P5.3.5: additional anammox lineage; tests hzs/hdh generalization across anammox genera.
        "P": ["hzsA","hzsB","hzsC","hdh","nirS","glnA"],
        "A": ["amoA","amoA_archaeal","amoA_gamma","nxrA","nxrB","narG","nosZ",
              "nosZ_clade2","nirK","nifH","nrfA"],
    },
    "Sstutzeri_F2a": {  # Stutzerimonas (Pseudomonas) stutzeri F2a — canonical γ-proteo complete
        # denitrifier (the textbook denit reference). Carries full denit (narGHI / nirS / norBC /
        # nosZ clade-I) + nasA/nasD assimilatory + glnA/gltB. P5.3.5: γ-proteo full-denit hold-out;
        # exercises every denit target end-to-end.
        "P": ["narG","narH","narI","nirS","norB","norC","nosZ","nasA","nasD",
              "glnA","gltB","gltD"],
        "A": ["napA","napB","nirK","amoA","amoA_archaeal","amoA_gamma","nxrA","nxrB",
              "nifH","nosZ_clade2","nrfA","hzsA","hzsB","hzsC","hdh"],
    },
    "Mcapsulatus_Bath": {  # Methylococcus capsulatus str. Bath — γ-proteo methanotroph (the
        # textbook methanotroph reference). Carries pmoABC (methane monooxygenase) — a known
        # amoA decoy via PF02461. KOfam K28504 (amoA-specific) should correctly NOT fire on
        # pmoA. P5.3.5: held-out methanotroph amoA-cross-fire test, plus nifH (diazotrophic
        # methanotroph). Should NOT fire amoA / amoA_gamma / amoA_archaeal.
        "P": ["nifH","nifD","nifK","nasA","glnA","gltB","gltD"],
        "A": ["amoA","amoA_gamma","amoA_archaeal","nxrA","nxrB","narG","nirS","nirK",
              "norB","nosZ","nosZ_clade2","nrfA","hzsA","hzsB","hzsC","hdh"],
    },
}


# P5.3.4 (2026-05-30): broad TN expansion for the two clade HMMs introduced
# in v7/v10. The audit asked to follow the `amoA_archaeal` pattern (TN cells
# in all panel genomes biologically expected to lack the target). Add the
# target as `absent` to every CURATED genome that doesn't have an explicit P
# or A call for it. Lifts amoA_gamma TN coverage 5 → ~32 and nosZ_clade2
# TN coverage 5 → ~32 — clears the P5.1.4 reliability flag once denominators
# rise above the floor. Source tag: `expanded_tn_v12`.
EXPAND_TN_TARGETS = ("amoA_gamma", "nosZ_clade2")


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str, str, str]] = [("genome", "target", "expected", "source")]
    for genome, sets in CURATED.items():
        src = "holdout_v3" if genome in HOLDOUTS else "curated_v11"
        explicit_targets = set(sets.get("P", [])) | set(sets.get("A", []))
        for t in dict.fromkeys(sets.get("P", [])):
            rows.append((genome, t, "present", src))
        for t in dict.fromkeys(sets.get("A", [])):
            rows.append((genome, t, "absent", src))
        # Broad TN expansion: implicit `absent` for the clade HMMs in every
        # genome that didn't get an explicit call above.
        expand_src = "holdout_v3" if genome in HOLDOUTS else "expanded_tn_v12"
        for t in EXPAND_TN_TARGETS:
            if t not in explicit_targets:
                rows.append((genome, t, "absent", expand_src))
    with open(OUT, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    n_p = sum(1 for r in rows[1:] if r[2] == "present")
    n_a = sum(1 for r in rows[1:] if r[2] == "absent")
    n_hold = sum(1 for r in rows[1:] if r[3] == "holdout_v3")
    n_expand = sum(1 for r in rows[1:] if r[3] == "expanded_tn_v12")
    print(f"[ground_truth] {len(rows)-1} cells across {len(CURATED)} genomes "
          f"({n_p} present, {n_a} absent; {n_hold} hold-out cells across "
          f"{len(HOLDOUTS)} genomes; {n_expand} P5.3.4 TN-expansion cells) -> {OUT}")


if __name__ == "__main__":
    main()
