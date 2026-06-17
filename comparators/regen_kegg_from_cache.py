#!/usr/bin/env python3
"""
Regenerate METABOLIC's per-genome present-KO set directly from the cached
intermediate_files/Hmmsearch_Outputs/*.hmmsearch_result.txt tblout files.

Rationale: METABOLIC ran hmmsearch with -T/--domT (thresholds applied at search
time), so each tblout already contains only passing protein hits, one line per
qualifying sequence. KEGG_identifier_result is just a per-genome reformat of
those counts. The expensive KEGG-module-completeness stage that runs *before*
that file is written is irrelevant to the (genome, KO) benchmark signal, so we
skip it and reconstruct the KO calls from the cached search output.

KO derivation per tblout: filename KO if it matches K#####, else the query-name
column (col3) if it matches K##### (covers custom HMMs like amoA->K10944).
Non-KO custom markers (aclA_alignment, aioA, ...) are dropped — they never
appear as K-prefixed lines in KEGG_identifier_result and the benchmark ignores
them.

--validate: compare the reconstructed present-KO sets against an existing
KEGG_identifier_result dir (the pilot's known-good output) and report agreement.
"""
import csv, glob, re, sys
from pathlib import Path
from collections import defaultdict

KO_RE = re.compile(r"^(K\d{5})")


def build_seqid2genome(proteome_dir: Path) -> dict:
    m = {}
    for fp in sorted(glob.glob(str(proteome_dir / "*.faa"))):
        gn = Path(fp).name[:-4]  # strip .faa
        with open(fp) as fh:
            for line in fh:
                if line.startswith(">"):
                    sid = line[1:].split()[0]
                    m[sid] = gn
    return m


def ko_for_tblout(path: Path) -> str:
    """KO from the filename if K#####, else from the query-name column."""
    base = path.name.split(".hmm.")[0]
    mm = KO_RE.match(base)
    if mm:
        return mm.group(1)
    # custom-named HMM: read the query name (col3) from the first hit line
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) >= 3:
                qm = KO_RE.match(cols[2])
                return qm.group(1) if qm else ""
    return ""


def regen(hmm_dir: Path, seqid2genome: dict):
    # present[(genome, ko)] = hit count
    present = defaultdict(int)
    genomes = set(seqid2genome.values())
    files = sorted(glob.glob(str(hmm_dir / "*.hmmsearch_result.txt")))
    for fp in files:
        p = Path(fp)
        ko = ko_for_tblout(p)
        if not ko:
            continue
        with open(fp) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                target = line.split()[0]
                gn = seqid2genome.get(target)
                if gn is None:
                    continue
                present[(gn, ko)] += 1
    return present, genomes, len(files)


def load_existing_kegg(kegg_dir: Path) -> dict:
    """Read an existing KEGG_identifier_result dir -> set of present (genome, ko)."""
    s = set()
    for fp in glob.glob(str(kegg_dir / "*.result.txt")):
        gn = Path(fp).name[: -len(".result.txt")]
        with open(fp) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0].startswith("K") and parts[1].strip():
                    s.add((gn, parts[0]))
    return s


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm-dir", type=Path, required=True,
                    help="metabolic_out/intermediate_files/Hmmsearch_Outputs")
    ap.add_argument("--proteome-dir", type=Path, required=True,
                    help="metabolic_in (the *.faa inputs)")
    ap.add_argument("--out", type=Path, help="output normalized (genome, ko) TSV of present calls")
    ap.add_argument("--validate", type=Path,
                    help="existing KEGG_identifier_result dir to compare against")
    args = ap.parse_args()

    seqid2genome = build_seqid2genome(args.proteome_dir)
    print(f"[regen] {len(set(seqid2genome.values()))} genomes, {len(seqid2genome)} proteins mapped",
          file=sys.stderr)
    present, genomes, nfiles = regen(args.hmm_dir, seqid2genome)
    pset = set(present.keys())
    print(f"[regen] scanned {nfiles} tblouts -> {len(pset)} present (genome,KO) cells",
          file=sys.stderr)

    if args.validate:
        gold = load_existing_kegg(args.validate)
        # restrict to genomes both sides cover
        gold_g = {g for g, _ in gold}
        mine_g = {g for g, _ in pset}
        common = gold_g & mine_g
        g = {x for x in gold if x[0] in common}
        m = {x for x in pset if x[0] in common}
        inter = g & m
        only_gold = g - m
        only_mine = m - g
        union = g | m
        jac = len(inter) / len(union) if union else 1.0
        print(f"\n[validate] genomes compared: {len(common)}")
        print(f"[validate] gold present cells:  {len(g)}")
        print(f"[validate] regen present cells: {len(m)}")
        print(f"[validate] agree (intersection): {len(inter)}")
        print(f"[validate] only in gold (regen MISSED): {len(only_gold)}")
        print(f"[validate] only in regen (regen EXTRA): {len(only_mine)}")
        print(f"[validate] Jaccard: {jac:.4f}")
        # which KOs drive disagreement
        from collections import Counter
        miss_ko = Counter(k for _, k in only_gold)
        extra_ko = Counter(k for _, k in only_mine)
        print(f"[validate] top MISSED KOs: {miss_ko.most_common(12)}")
        print(f"[validate] top EXTRA KOs:  {extra_ko.most_common(12)}")

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["genome", "ko"], delimiter="\t")
            w.writeheader()
            for (gn, ko) in sorted(pset):
                w.writerow({"genome": gn, "ko": ko})
        print(f"[regen] wrote {args.out} ({len(pset)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
