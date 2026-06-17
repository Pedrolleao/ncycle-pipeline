#!/usr/bin/env python3
"""Supplemental NxrA: explicit curated NOB accessions. Many NOB NXR alpha subunits
are UniProt-annotated as 'nitrate reductase (quinone)' (the exact NxrA/NarG ambiguity
this reference set exists to resolve). We fetch them by accession and label as nxr,
recording the annotation caveat per entry.
"""
import json, sys, time, urllib.request, os
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# acc, subtype, organism_label, note
CURATED = [
    # Nitrobacter NXR alpha (annotated 'nitrate reductase (quinone)' in UniProt;
    # Nitrobacter is an obligate NOB whose Mo-enzyme is NXR, Starkenburg 2006 Nwi genome)
    ("Q3SUK2", "Nitrobacter-type", "Nitrobacter winogradskyi Nb-255", "UniProt-annot='nitrate reductase(quinone)'; NOB NXR alpha"),
    ("Q1QPP1", "Nitrobacter-type", "Nitrobacter hamburgensis X14", "UniProt-annot='nitrate reductase(quinone)'; NOB NXR alpha"),
    ("A0ACC6ADH6", "Nitrobacter-type", "Nitrobacter winogradskyi", "UniProt-annot='nitrate reductase(quinone)'; NOB NXR alpha"),
    # Nitrococcus mobilis NXR (Gammaproteobacteria NOB)
    ("A4BM17", "Nitrococcus-type", "Nitrococcus mobilis Nb-231", "UniProt-annot='nitrate reductase(quinone)'; NOB NXR alpha"),
    # Nitrospina (deep-branching marine NOB) - explicitly nxrA-annotated
    ("A0ABM9HBT2", "Nitrospina-type", "Nitrospina watsonii", "nxrA-annotated; periplasmic/divergent marine NOB NXR"),
    # Nitrospira defluvii NXR alpha (periplasmic NXR; Lucker 2010) - try canonical
    ("D8PD96", "Nitrospira-type", "Candidatus Nitrospira defluvii", "periplasmic Nitrospira-type NXR alpha (Lucker 2010)"),
    # Nitrospira moscoviensis
    ("A0A0F2QU07", "Nitrospira-type", "Nitrospira moscoviensis", "Nitrospira-type NXR alpha"),
    # comammox Nitrospira inopinata
    ("A0A1H6QXC0", "Nitrospira-type", "Candidatus Nitrospira inopinata", "comammox NXR alpha"),
    # Nitrolancea (already have one; add a 2nd Chloroflexi-type if distinct) - skip dup
    # Kuenenia anammox NXR (NarGH-like, runs nitrite oxidation in anammox)
    ("Q1Q0V8", "anammox-type", "Candidatus Kuenenia stuttgartiensis", "anammox NXR/Nar-like alpha (nitrite oxidation to nitrate)"),
]

def fetch_json(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    req = urllib.request.Request(url, headers={"User-Agent":"ncycle-curation/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=40) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        sys.stderr.write(f"FAIL {acc}: {e}\n"); return None

def main():
    rows = []
    for acc, subtype, org, note in CURATED:
        d = fetch_json(acc); time.sleep(0.4)
        if not d:
            print("DROP", acc); continue
        seq = d.get("sequence",{}).get("value","")
        L = d.get("sequence",{}).get("length",0)
        rev = "reviewed" in d.get("entryType","").lower()
        if not seq or not (900 <= L <= 1350):
            print("DROP-len", acc, L); continue
        realorg = d.get("organism",{}).get("scientificName", org)
        rows.append((acc, "nxrA", "nxr", subtype, org, "UniProt",
                     f"len={L};{'reviewed' if rev else 'unreviewed'};{note}", seq))
        print("OK", acc, L, subtype)
    json.dump(rows, open(os.path.join(OUTDIR,"_nxr_supp.json"),"w"))

if __name__ == "__main__":
    main()
