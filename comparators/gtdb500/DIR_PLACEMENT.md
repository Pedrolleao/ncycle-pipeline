# ncycle — nxrA/narG identity by PLACEMENT + OPERON SYNTENY (B9.4)

Two orthogonal signals best-hit (B9.3) lacked: a combined ref+query ML tree (IQ-TREE MFP + 1000 UFBoot, midpoint-rooted; clade-membership classification) and operon synteny (narGHJI vs nxr gene neighborhood from Prodigal ordinals). Synteny is decisive where present; placement otherwise.


## Method self-validation (characterized genomes)

Verdict matches ncycle on **29/30** characterized calls (97%).


## Adjudication

- nxrA/narG calls with a definite verdict: **54**; **51** (94%) agree with ncycle.
- Calls resolved by **operon synteny**: **41**.
- **Best-hit disagreements now adjudicated:** of the 15 calls where best-hit contradicted ncycle, **12** are confirmed in ncycle's favour by placement+synteny — i.e. best-hit was mislabeling the Nitrobacter-type NXR≈NarG paralogy, as predicted; ncycle's HMM+gating call stands.


## Likely ncycle over-calls flagged by BOTH orthogonal signals

**3** calls where placement AND best-hit agree against ncycle (candidate false positives for curation):
- `GCF_003176655_1` (denit_Pseudomonas): ncycle **nxrA**, but placement=nar + best-hit=nar ⇒ **nar**.
- `GCA_031960865_1` (p__Campylobacterota_A): ncycle **nxrA**, but placement=nar + best-hit=nar ⇒ **nar**.
- `GCA_050253135_1` (p__Campylobacterota_A): ncycle **nxrA**, but placement=nar + best-hit=nar ⇒ **nar**.

The orthogonal check works both ways — it vindicates ncycle on the Nitrobacter paralogy AND surfaces genuine over-calls (notably an nxrA call on a *Pseudomonas* denitrifier, nar by both signals) for follow-up.


## Per-call detail

Full table: `DIR_PLACEMENT.tsv`.

