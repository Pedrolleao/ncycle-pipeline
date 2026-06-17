#!/usr/bin/env python3
"""Build the combined nxr/nar ref+query FASTA (one tree) with safe leaf names + label map."""
import csv, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GT = HERE.parents[1] / "comparators" / "gtdb500"
PMAP = {r["sample"]: r for r in csv.DictReader(open(GT / "panel_map.tsv"), delimiter="\t")}
NCY = {(r["genome"], r["ncycle_gene"]): r for r in csv.DictReader(open(GT / "DIR_ACCURACY.tsv"), delimiter="\t")}


def read_faa(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


refs = read_faa(HERE / "nxr_nar_refs.faa")
qrys = read_faa(GT / "_nxrnar_query.faa")
label, out = {}, HERE / "trees" / "nxrnar_in.faa"
out.parent.mkdir(exist_ok=True)
ri = qi = 0
with open(out, "w") as fh:
    for name, seq in refs.items():
        acc, gene, typ, org = name.split("|")
        lid = f"R{ri:03d}"; ri += 1
        label[lid] = {"kind": "ref", "type": typ, "gene": gene, "acc": acc, "org": org}
        fh.write(f">{lid}\n{seq}\n")
    for name, seq in qrys.items():
        genome, tgt, pid = name.split("__")
        lid = f"Q{qi:03d}"; qi += 1
        cl = PMAP.get(genome, {}).get("clade", "")
        label[lid] = {"kind": "query", "genome": genome, "ncycle_gene": tgt,
                      "ncycle_type": "nxr" if tgt == "nxrA" else "nar",
                      "clade": cl, "characterized": (not cl.startswith("p__")),
                      "protein_id": pid}
        fh.write(f">{lid}\n{seq}\n")
json.dump(label, open(HERE / "trees" / "label_map.json", "w"), indent=0)
print(f"nxr/nar: {ri} refs + {qi} queries -> {out}", file=sys.stderr)
