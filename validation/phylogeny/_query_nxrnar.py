#!/usr/bin/env python3
"""Query UniProt for NxrA (nitrite oxidoreductase alpha) and NarG (membrane-bound
nitrate reductase alpha). Type-II DMSO reductase Mo-bis-PGD superfamily paralogs.
NxrA ~1100-1300 aa; NarG ~1200-1250 aa. Strict name+length filter, prefer reviewed.
"""
import json, sys, time, urllib.request, urllib.parse, os
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# (organism, gene nxrA|narG, type nxr|nar, subtype, note)
TARGETS = [
    # ---- NxrA (nitrite oxidizers) ----
    ("Nitrobacter winogradskyi", "nxrA", "nxr", "Nitrobacter-type"),
    ("Nitrobacter hamburgensis", "nxrA", "nxr", "Nitrobacter-type"),
    ("Nitrococcus mobilis", "nxrA", "nxr", "Nitrococcus-type"),
    ("Candidatus Nitrospira defluvii", "nxrA", "nxr", "Nitrospira-type"),
    ("Nitrospira", "nxrA", "nxr", "Nitrospira-type"),
    ("Nitrolancea hollandica", "nxrA", "nxr", "Nitrolancea-type"),
    ("Nitrospina gracilis", "nxrA", "nxr", "Nitrospina-type"),
    ("Candidatus Nitrotoga", "nxrA", "nxr", "Nitrotoga-type"),
    ("Candidatus Kuenenia stuttgartiensis", "nxrA", "nxr", "anammox-type"),
    ("Candidatus Brocadia", "nxrA", "nxr", "anammox-type"),
    # ---- NarG (denitrifiers / nitrate reducers) ----
    ("Escherichia coli", "narG", "nar", "canonical_NarG"),
    ("Paracoccus denitrificans", "narG", "nar", "canonical_NarG"),
    ("Pseudomonas aeruginosa", "narG", "nar", "canonical_NarG"),
    ("Bacillus subtilis", "narG", "nar", "canonical_NarG"),
    ("Thermus thermophilus", "narG", "nar", "canonical_NarG"),
    ("Mycobacterium tuberculosis", "narG", "nar", "canonical_NarG"),
    ("Rhodobacter sphaeroides", "narG", "nar", "canonical_NarG"),
    ("Halomonas", "narG", "nar", "canonical_NarG"),
    ("Salmonella", "narG", "nar", "canonical_NarG"),
    ("Geobacillus", "narG", "nar", "canonical_NarG"),
    ("Corynebacterium glutamicum", "narG", "nar", "canonical_NarG"),
    ("Magnetospirillum", "narG", "nar", "canonical_NarG"),
]

def uquery(q):
    fields = "accession,gene_names,protein_name,organism_name,length,reviewed,sequence"
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urllib.parse.urlencode({"query": q, "fields": fields, "format": "json", "size": "40"}))
    req = urllib.request.Request(url, headers={"User-Agent":"ncycle-curation/1.0"})
    for a in range(3):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            sys.stderr.write(f"retry {a}: {e}\n"); time.sleep(2)
    return {"results": []}

def candidates(res, gene):
    out = []
    for r in res.get("results", []):
        L = r.get("sequence",{}).get("length",0)
        seq = r.get("sequence",{}).get("value","")
        pd = r.get("proteinDescription",{})
        pname = (pd.get("recommendedName",{}).get("fullName",{}).get("value","")) or ""
        if not pname:
            for s in pd.get("submissionNames",[]):
                pname = pname or s.get("fullName",{}).get("value","")
        gn = " ".join(g.get("geneName",{}).get("value","") for g in r.get("genes",[])).lower()
        acc = r.get("primaryAccession","")
        rev = "reviewed" in r.get("entryType","").lower()
        nm = pname.lower()
        org = r.get("organism",{}).get("scientificName","")
        if not seq or not (900 <= L <= 1350):
            continue
        if gene == "nxrA":
            is_match = ("nitrite oxidoreductase" in nm or "nxr" in gn or
                        ("nitrite" in nm and "oxidoreductase" in nm))
            alpha = "alpha" in nm or "nxra" in gn or gn.strip()=="nxra"
        else:
            is_match = (("nitrate reductase" in nm and "alpha" in nm) or "narg" in gn or
                        ("respiratory nitrate reductase" in nm))
            alpha = "alpha" in nm or "narg" in gn
            # exclude periplasmic NapA and assimilatory NasA
            if "periplasmic" in nm or "nap" in gn or "assimilatory" in nm or "nas" in gn:
                continue
        if not is_match:
            continue
        score = (10 if rev else 0) + (5 if alpha else 0)
        out.append((score, acc, org, L, seq, pname or gn, rev))
    out.sort(reverse=True)
    return out

def main():
    rows, seen, log = [], set(), []
    for orgq, gene, typ, subtype in TARGETS:
        gfilt = "nxrA" if gene == "nxrA" else "narG"
        if gene == "nxrA":
            q = f'(organism_name:"{orgq}") AND (gene:nxrA OR protein_name:"nitrite oxidoreductase")'
        else:
            q = f'(organism_name:"{orgq}") AND (gene:narG OR protein_name:"nitrate reductase alpha" OR protein_name:"respiratory nitrate reductase")'
        res = uquery(q); time.sleep(0.4)
        cands = candidates(res, gene)
        picked = next((c for c in cands if c[1] not in seen), None)
        if picked is None:
            log.append(f"NONE {orgq} {gene}"); continue
        sc, acc, org, L, seq, pname, rev = picked
        seen.add(acc)
        rows.append((acc, gene, typ, subtype, org, "UniProt",
                     f"len={L};{'reviewed' if rev else 'unreviewed'};{pname[:45]}", seq))
        log.append(f"OK {orgq} {gene} -> {acc} L={L} rev={rev}")
    json.dump(rows, open(os.path.join(OUTDIR,"_nxrnar_seqs.json"),"w"))
    open(os.path.join(OUTDIR,"_nxrnar_log.txt"),"w").write("\n".join(log))
    from collections import Counter
    print("rows", len(rows), dict(Counter((g,t) for _,g,t,_,_,_,_,_ in rows)))
    print("\n".join(log))

if __name__ == "__main__":
    main()
