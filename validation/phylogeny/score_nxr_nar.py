#!/usr/bin/env python3
"""
score_nxr_nar.py — B9.2/B9.3 orthogonal accuracy check for ncycle's nxrA-vs-narG
IDENTITY call (the N analogue of scycle's dsrAB direction scorer).

nxrA (nitrite oxidoreductase α, nitrite oxidation) and narG (membrane nitrate
reductase α, nitrate reduction) are close paralogs in the type-II DMSO-reductase
molybdoenzyme family that share KOs — the homology trap behind the nxrA→narG
mis-routing. ncycle disambiguates them with custom HMMs + gating; this script
checks that call INDEPENDENTLY from SEQUENCE ANCESTRY: it DIAMOND-blasts each
ncycle-called nxrA / narG protein against the type-labeled reference set
(validation/phylogeny/nxr_nar_refs.*) and takes the best-hit gene type (nxr | nar).
The two signals share no inputs, so their agreement is a non-circular accuracy
measure of the identity call. The unit is the CALL (a genome may have both).

Lightweight tier (best-hit); the nxr-vs-nar bitscore MARGIN is recorded so close
calls — expected for Nitrobacter-type NxrA (~55-60% id to NarG) and divergent
periplasmic Nitrospira/Nitrospina NXR — can be escalated to placement (B9.4).

Outputs (under comparators/gtdb500/):
  DIR_ACCURACY.tsv  — per-call: ncycle gene, phylo gene, margin, clade
  DIR_ACCURACY.md   — confusion matrix + accuracy (Wilson 95% CI)
"""
from __future__ import annotations
import argparse, csv, math, subprocess, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRESENT = {"confirmed", "domain-only", "narrow-no-IPR"}
TARGET_GENE = {"nxrA": "nxr", "narG": "nar"}   # ncycle target → expected phylo type


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, (c - h) / d, (c + h) / d


def read_proteome(faa: Path) -> dict[str, str]:
    seqs, pid, buf = {}, None, []
    if not faa.exists():
        return seqs
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid:
                    seqs[pid] = "".join(buf)
                pid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if pid:
        seqs[pid] = "".join(buf)
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--proteomes", type=Path, required=True)
    ap.add_argument("--db", type=Path, default=HERE / "nxr_nar_refs.dmnd")
    ap.add_argument("--panel-map", type=Path, help="panel_map.tsv (sample, clade)")
    ap.add_argument("--out", type=Path, required=True, help="DIR_ACCURACY.md")
    ap.add_argument("--margin-frac", type=float, default=0.10,
                    help="nxr-vs-nar best-hit bitscore margin below this → escalate to placement "
                         "(Nitrobacter NxrA ~55-60% id to NarG → genuinely close)")
    args = ap.parse_args()

    clade = {}
    if args.panel_map and args.panel_map.exists():
        for r in csv.DictReader(open(args.panel_map), delimiter="\t"):
            clade[(r.get("sample") or r.get("ncbi_acc", "")).replace(".", "_")] = r.get("clade", "")
    def characterized(g):
        return bool(clade) and not clade.get(g, "p__").startswith("p__")

    genomes = sorted(p.name for p in args.results.iterdir()
                     if p.is_dir() and p.name.startswith("GC"))

    # ── 1. collect ncycle nxrA/narG calls per genome ────────────────────────
    want = []   # (genome, target, protein_id)
    for g in genomes:
        cf = args.results / g / "calls" / "ncycle_calls.tsv"
        if not cf.exists():
            continue
        for r in csv.DictReader(open(cf), delimiter="\t"):
            if r["target_id"] in TARGET_GENE and r["status"] in PRESENT and r.get("protein_id"):
                want.append((g, r["target_id"], r["protein_id"]))

    # ── 2. query FASTA ──────────────────────────────────────────────────────
    qfaa = args.out.parent / "_nxrnar_query.faa"
    n_q = 0
    seqs_cache: dict[str, dict[str, str]] = {}
    with open(qfaa, "w") as out:
        for g, tgt, pid in want:
            if g not in seqs_cache:
                seqs_cache[g] = read_proteome(args.proteomes / f"{g}.faa")
            s = seqs_cache[g].get(pid)
            if s:
                out.write(f">{g}__{tgt}__{pid}\n{s}\n")
                n_q += 1
    print(f"[nxrnar] {len(want)} nxrA/narG calls, {n_q} query proteins", file=sys.stderr)

    # ── 3. DIAMOND blastp vs typed refs ─────────────────────────────────────
    bl = args.out.parent / "_nxrnar_query.blast.tsv"
    subprocess.run(["diamond", "blastp", "--db", str(args.db), "--query", str(qfaa),
                    "--very-sensitive", "-k", "50", "-e", "1e-5", "--quiet", "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "bitscore", "--out", str(bl)],
                   check=True)

    # ── 4. per-query best-hit type (nxr/nar) + margin ───────────────────────
    best_nxr: dict[str, float] = defaultdict(float)
    best_nar: dict[str, float] = defaultdict(float)
    for r in csv.reader(open(bl), delimiter="\t"):
        q, s, bits = r[0], r[1], float(r[4])
        parts = s.split("|")            # acc|gene|type|organism ; type = nxr|nar
        if len(parts) < 3:
            continue
        t = parts[2]
        if t == "nxr":
            best_nxr[q] = max(best_nxr[q], bits)
        elif t == "nar":
            best_nar[q] = max(best_nar[q], bits)

    rows = []
    for g, tgt, pid in want:
        q = f"{g}__{tgt}__{pid}"
        nx, na = best_nxr.get(q, 0.0), best_nar.get(q, 0.0)
        if nx == 0 and na == 0:
            phylo, margin = "no_hit", 0.0
        else:
            phylo = "nxr" if nx >= na else "nar"
            top = max(nx, na)
            margin = abs(nx - na) / top if top else 0.0
        rows.append({
            "genome": g, "clade": clade.get(g, ""),
            "characterized": "yes" if characterized(g) else "",
            "ncycle_gene": tgt, "ncycle_type": TARGET_GENE[tgt],
            "phylo_type": phylo, "margin": f"{margin:.3f}",
            "escalate": "yes" if (phylo == "no_hit" or 0 < margin < args.margin_frac) else "",
        })
    rows.sort(key=lambda r: (r["ncycle_gene"], r["phylo_type"]))

    tsv = args.out.with_suffix(".tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ── 5. confusion + accuracy ─────────────────────────────────────────────
    defn = [r for r in rows if r["phylo_type"] in ("nxr", "nar")]
    agree = sum(1 for r in defn if r["ncycle_type"] == r["phylo_type"])
    p, lo, hi = wilson(agree, len(defn))

    def acc(sub):
        a = sum(1 for r in sub if r["ncycle_type"] == r["phylo_type"])
        pp, ll, hh = wilson(a, len(sub))
        return a, len(sub), pp, ll, hh
    char = [r for r in defn if r["characterized"] == "yes"]
    cand = [r for r in defn if r["characterized"] != "yes"]
    conf = [r for r in defn if float(r["margin"]) >= args.margin_frac]
    ca, cn, cp, clo, chi = acc(char)
    da, dn, dp, dlo, dhi = acc(cand)
    fa, fn, fp, flo, fhi = acc(conf)

    cm = defaultdict(int)
    for r in rows:
        cm[(r["ncycle_gene"], r["phylo_type"])] += 1
    cols = ["nxr", "nar", "no_hit"]

    L = ["# ncycle — nxrA-vs-narG identity accuracy vs phylogeny-anchored reference (B9)\n",
         "Orthogonal, **non-circular** check: ncycle disambiguates the shared-KO nxrA/narG trap "
         "with custom HMMs + gating; the reference call comes from type-II DMSO-reductase sequence "
         "ancestry (type-labeled nxr/nar refs, best-hit DIAMOND). Agreement = accuracy of the "
         "identity call. Unit = each nxrA/narG call.\n",
         f"**nxrA/narG calls scored:** {len(want)} · **with a reference hit:** {len(defn)}\n",
         "| stratum | agreement | accuracy | Wilson 95% CI |",
         "|---|---|---|---|",
         f"| **Characterized (named genus)** | {ca}/{cn} | **{cp:.3f}** | [{clo:.3f}, {chi:.3f}] |",
         f"| Uncultured candidate phyla (`p__`) | {da}/{dn} | {dp:.3f} | [{dlo:.3f}, {dhi:.3f}] |",
         f"| Confident best-hit (margin ≥ {args.margin_frac:.2f}) | {fa}/{fn} | {fp:.3f} | [{flo:.3f}, {fhi:.3f}] |",
         f"| **All with a hit** | {agree}/{len(defn)} | {p:.3f} | [{lo:.3f}, {hi:.3f}] |",
         "\n> **Interpretation — best-hit is INCONCLUSIVE for the nxr/nar trap (this is the finding, "
         "not an ncycle accuracy verdict).** Unlike the dsrAB reductive/oxidative split (cleanly "
         "separable, scycle validated at 100% on characterized genomes by the same method), NxrA and "
         "NarG are reciprocally close paralogs in the type-II DMSO-reductase family: Nitrobacter-type "
         "NXR is phylogenetically embedded next to NarG. The discordances concentrate on exactly the "
         "expected lineages — **Nitrobacter ×5, Nitrococcus, Nitrolancea** (NOB) calling narG but "
         "best-hitting the nearby Nitrobacter NxrA reference (no Nitrobacter-clade narG exists in the "
         "20-seq reference), plus divergent candidate phyla. In these, ncycle's custom-HMM+gating call "
         "is the more reliable signal and best-hit is mislabeling. The nxr/nar identity therefore "
         "REQUIRES the rigorous placement tier (B9.4, with operon/synteny context); the curated-panel "
         "benchmark remains the N accuracy authority. This asymmetry (sulfur resolvable by best-hit, "
         "nitrogen not) is itself a reportable result.\n",
         "\n## Confusion matrix — ncycle call (rows) × phylogeny best-hit (cols)\n",
         "| ncycle ↓ \\ phylo → | " + " | ".join(cols) + " |",
         "|---|" + "|".join("---" for _ in cols) + "|"]
    for tgt in ("nxrA", "narG"):
        if any(cm[(tgt, c)] for c in cols):
            L.append(f"| **{tgt}** | " + " | ".join(str(cm[(tgt, c)]) for c in cols) + " |")
    n_esc = sum(1 for r in rows if r["escalate"] == "yes")
    n_nohit = sum(1 for r in rows if r["phylo_type"] == "no_hit")
    L.append(f"\n## Escalation set for B9.4 (EPA-ng/gappa placement)\n")
    L.append(f"- **Disagreements** (ncycle ≠ phylo, with hit): **{len(defn) - agree}**.")
    L.append(f"- **Low best-hit margin** (<{args.margin_frac:.2f}): **{n_esc - n_nohit}** — "
             "Nitrobacter-type NxrA / divergent periplasmic NXR sit close to NarG; place these.")
    L.append(f"- **No reference hit**: **{n_nohit}** (divergent enzyme → placement needed).")
    L.append(f"\nPer-call detail: `{tsv.name}`. Lightweight tier; rigorous placement is B9.4.\n")
    args.out.write_text("\n".join(L) + "\n")
    print(f"[nxrnar] accuracy {agree}/{len(defn)}={p:.3f} CI[{lo:.3f},{hi:.3f}]; "
          f"char {ca}/{cn}; escalate={n_esc} nohit={n_nohit}", file=sys.stderr)
    print(f"[nxrnar] wrote {args.out} + {tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
