#!/usr/bin/env python3
import json, os
from collections import Counter
OUTDIR = os.path.dirname(os.path.abspath(__file__))
main = json.load(open(os.path.join(OUTDIR,"_nxrnar_seqs.json")))
supp = json.load(open(os.path.join(OUTDIR,"_nxr_supp.json")))

# ncycle presence-gate refs (circularity overlap) from config/targets.yaml
NCYCLE_PRESENCE = {"Q71RT9","A0A0S4LQF4","A0A7S8FGD0",  # nxrA gate
                   "Q51075","A0A0S4LK29","M1KVL1"}      # nxrB gate (different gene but note)

rows = main + supp
seen = set()
faa, tsv, kept = [], ["accession\tgene\ttype\tsubtype\torganism\tsource\tnotes"], []
for acc, gene, typ, subtype, org, src, note, seq in rows:
    if acc in seen:
        continue
    seen.add(acc)
    seq = seq.replace("*","").strip()
    if len(seq) < 300:
        continue
    overlap = " [ncycle-presence-ref]" if acc in NCYCLE_PRESENCE else ""
    org_clean = org.split(" (strain")[0]
    faa.append(f">{acc}|{gene}|{typ}|{org_clean.replace(' ','_')}")
    for i in range(0,len(seq),60):
        faa.append(seq[i:i+60])
    tsv.append(f"{acc}\t{gene}\t{typ}\t{subtype}\t{org_clean}\t{src}\t{note}{overlap}")
    kept.append((gene,typ,subtype,acc in NCYCLE_PRESENCE))

open(os.path.join(OUTDIR,"nxr_nar_refs.faa"),"w").write("\n".join(faa)+"\n")
open(os.path.join(OUTDIR,"nxr_nar_refs.tsv"),"w").write("\n".join(tsv)+"\n")
print("total", len(kept))
print("gene x type", dict(Counter((g,t) for g,t,_,_ in kept)))
print("subtype", dict(Counter(s for _,_,s,_ in kept)))
print("ncycle-overlap", sum(1 for *_,o in kept if o))

if __name__ == "__main__":
    pass
