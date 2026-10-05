"""The nitrogen-cycle diagram: compounds, reaction steps, and how a genome's
calls light the steps up.

Single source of truth for the static cycle maps (make_cycle_map.py) and
the interactive one in report.html (make_html_report.py ships the geometry
computed here as JSON, so the browser only draws polylines).

A step is `complete` when every gene of at least one route is present
(confirmed or domain-only — the same rule as complex completeness), `partial`
when some route gene is present but no route is whole, otherwise `absent`.
This is a drawing-level summary; the authoritative per-process calls remain
complex_completeness.tsv and synergy_completeness.tsv.
"""

from __future__ import annotations

import math

# Canvas (data units, y up).
XLIM = (-17.0, 100.0)
YLIM = (-5.0, 97.0)

# id: (label, x, y, width, height)
NODES = {
    "N2":    ("N₂",        52.0, 88.0,  9.0, 7.0),
    "NH4":   ("NH₄⁺",      20.0, 56.0, 11.0, 7.0),
    "NH2OH": ("NH₂OH",     29.0, 30.0, 15.0, 7.0),
    "NO2":   ("NO₂⁻",      52.0, 23.0, 11.0, 7.0),
    "NO3":   ("NO₃⁻",      52.0,  1.5, 14.0, 7.0),
    "NO":    ("NO",        76.0, 33.0,  8.0, 7.0),
    "N2O":   ("N₂O",       84.0, 62.0, 10.0, 7.0),
    "N2H4":  ("N₂H₄",      51.0, 57.0, 12.0, 7.0),
    "ORG":   ("organic N", -5.0, 67.0, 19.0, 7.0),
    "UREA":  ("urea",      -5.0, 46.0, 11.0, 7.0),
}

# Short display names for the clade-specific targets.
DISPLAY = {
    "amoA_archaeal": "amoA-arch",
    "amoA_gamma":    "amoA-γ",
    "nosZ_clade2":   "nosZ-II",
}

_AMOA = ["amoA", "amoA_gamma", "amoA_archaeal"]

# arrows: (from, to, bow, (dx0, dy0), (dx1, dy1)) — bow is the offset of the
#   curve's control point to the LEFT of the direction of travel; the two
#   offsets shift the start / end anchor off the node centre (parallel arrows).
# routes: alternative gene sets, any one of which performs the step; a slot
#   that is a list is any-of (clade paralogs).
# lines:  how the genes are written next to the arrow (may include accessory
#   genes that no route requires, e.g. nrfH).
# label:  (x, y, ha) of the first label line; further lines stack downwards.
STEPS = [
    {"id": "fix", "name": "N₂ fixation", "pathway": "nitrogen_fixation",
     "arrows": [("N2", "NH4", -9.0, (0, 0), (0, 0))],
     "routes": [["nifH", "nifD", "nifK"]],
     "lines": [["nifH", "nifD", "nifK"]],
     "label": (26.5, 83.5, "right")},
    {"id": "amo", "name": "Ammonia oxidation", "pathway": "nitrification_ammonia",
     "arrows": [("NH4", "NH2OH", -4.0, (0, 0), (0, 0))],
     "routes": [[_AMOA, "amoB", "amoC"]],
     "lines": [[_AMOA, "amoB", "amoC"]],
     "label": (19.5, 38.5, "right")},
    {"id": "hao", "name": "Hydroxylamine oxidation", "pathway": "nitrification_ammonia",
     "arrows": [("NH2OH", "NO2", -3.0, (0, 0), (0, 0))],
     "routes": [["hao"]],
     "lines": [["hao"]],
     "label": (38.0, 20.0, "center")},
    {"id": "nxr", "name": "Nitrite oxidation", "pathway": "nitrification_nitrite",
     "arrows": [("NO2", "NO3", 0.0, (-4.0, 0), (-5.0, 0))],
     "routes": [["nxrA", "nxrB"]],
     "lines": [["nxrA", "nxrB"]],
     "label": (44.5, 12.5, "right")},
    {"id": "nar", "name": "Nitrate reduction (respiratory)", "pathway": "denitrification",
     "arrows": [("NO3", "NO2", 0.0, (0.5, 0), (0.5, 0))],
     "routes": [["narG", "narH", "narI"], ["napA", "napB"]],
     "lines": [["narG", "narH", "narI"], ["napA", "napB"]],
     "label": (60.0, 17.0, "left")},
    {"id": "nas", "name": "Nitrate reduction (assimilatory)", "pathway": "assimilatory",
     "arrows": [("NO3", "NO2", 0.0, (5.0, 0), (4.0, 0))],
     "routes": [["nasA"], ["narB"], ["NR"]],
     "lines": [["nasA", "nasB"], ["narB", "NR"]],
     "label": (60.0, 8.5, "left")},
    {"id": "nir", "name": "Nitrite reduction to NO", "pathway": "denitrification",
     "arrows": [("NO2", "NO", -3.0, (0, 0), (0, 0))],
     "routes": [["nirK"], ["nirS"]],
     "lines": [["nirK", "nirS"]],
     "label": (67.0, 24.0, "left")},
    {"id": "nor", "name": "NO reduction", "pathway": "denitrification",
     "arrows": [("NO", "N2O", -4.0, (0, 0), (0, 0))],
     "routes": [["norB", "norC"], ["norZ"]],
     "lines": [["norB", "norC"], ["norZ"]],
     "label": (86.0, 48.5, "left")},
    {"id": "nos", "name": "N₂O reduction", "pathway": "denitrification",
     "arrows": [("N2O", "N2", -9.0, (0, 0), (0, 0))],
     "routes": [["nosZ"], ["nosZ_clade2"]],
     "lines": [["nosZ", "nosZ_clade2"]],
     "label": (76.0, 83.5, "left")},
    {"id": "dnra", "name": "Nitrite reduction to ammonium (DNRA)", "pathway": "dnra",
     "arrows": [("NO2", "NH4", 4.5, (0.0, 1.5), (6.0, -1.0))],
     "routes": [["nrfA"], ["nirB", "nirD"]],
     "lines": [["nrfA", "nrfH"], ["nirB", "nirD"]],
     "label": (43.0, 47.0, "left")},
    {"id": "nia", "name": "Nitrite reduction to ammonium (assimilatory)", "pathway": "assimilatory",
     "arrows": [("NO2", "NH4", 1.0, (-3.5, -0.5), (-0.5, -2.0))],
     "routes": [["nirA"], ["nasD"]],
     "lines": [["nirA", "nasD"]],
     "label": (43.0, 38.0, "left")},
    {"id": "hzs", "name": "Hydrazine synthesis (anammox)", "pathway": "anammox",
     "arrows": [("NO", "N2H4", 0.0, (0, 0), (0, 0)),
                ("NH4", "N2H4", 0.0, (0, 0.5), (0, 0))],
     "routes": [["hzsA", "hzsB", "hzsC"]],
     "lines": [["hzsA", "hzsB", "hzsC"]],
     "label": (59.5, 56.5, "left")},
    {"id": "hdh", "name": "Hydrazine oxidation (anammox)", "pathway": "anammox",
     "arrows": [("N2H4", "N2", 0.0, (0, 0), (0, 0))],
     "routes": [["hdh"]],
     "lines": [["hdh"]],
     "label": (54.0, 72.5, "left")},
    {"id": "gs", "name": "Ammonia assimilation", "pathway": "ammonia_assimilation",
     "arrows": [("NH4", "ORG", 0.0, (-1.0, 1.5), (3.0, 0))],
     "routes": [["glnA", "gltB", "gltD"], ["gdhA"]],
     "lines": [["glnA", "gltB", "gltD"], ["gdhA"]],
     "label": (-5.0, 78.5, "center")},
    {"id": "ure", "name": "Urea hydrolysis", "pathway": "organic_n_mineralization",
     "arrows": [("UREA", "NH4", 0.0, (2.0, 0), (-1.0, -1.5))],
     "routes": [["ureA", "ureB", "ureC"]],
     "lines": [["ureA", "ureB", "ureC"]],
     "label": (-5.0, 53.0, "center")},
]

# Figure sizing used by make_cycle_map.py (inches / points).
FIG = {"single_w": 7.4, "grid_w": 5.3, "single_fs": 9.0, "grid_fs": 6.6,
       "grid_cols": 3}

LINE_STEP = 3.6          # distance between stacked label lines (data units)
HEAD_LEN, HEAD_W = 2.4, 1.25
NODE_PAD = 1.3           # gap between an arrow end and the node box


# ───────────────────────────── evaluation ────────────────────────────────────

def _slot_genes(slot) -> list[str]:
    return [slot] if isinstance(slot, str) else list(slot)


def _token(slot, codes: dict[str, int]) -> dict:
    """The gene to print for a slot: the best-supported alternative (first one
    when none is present)."""
    genes = _slot_genes(slot)
    rank = {2: 0, 1: 1, -1: 2, 0: 3}
    best = min(genes, key=lambda g: (rank[codes.get(g, 0)], genes.index(g)))
    return {"gene": best, "label": DISPLAY.get(best, best),
            "code": codes.get(best, 0)}


def genome_context(calls: dict[str, dict]) -> dict:
    """Genome-level facts that change how steps are drawn (e.g. the direction
    of a reversible enzyme). The nitrogen map has none."""
    return {}


def context_note(ctx: dict) -> str:
    """One-line description of the context, for subtitles ('' if none)."""
    return ""


def evaluate(codes: dict[str, int], ctx: dict | None = None) -> dict[str, dict]:
    """{step_id: {state, reversed, lines:[[token,…],…]}} for one genome."""
    out: dict[str, dict] = {}
    for st in STEPS:
        def present(slot) -> bool:
            return any(codes.get(g, 0) in (1, 2) for g in _slot_genes(slot))
        route_full = any(all(present(s) for s in r) for r in st["routes"])
        any_gene = any(present(s) for r in st["routes"] for s in r)
        state = "complete" if route_full else "partial" if any_gene else "absent"
        out[st["id"]] = {
            "state": state,
            "reversed": False,
            "lines": [[_token(s, codes) for s in line] for line in st["lines"]],
        }
    return out


# ───────────────────────────── geometry ──────────────────────────────────────

def _inside(pt, node, pad) -> bool:
    _, x, y, w, h = NODES[node]
    return abs(pt[0] - x) <= w / 2 + pad and abs(pt[1] - y) <= h / 2 + pad


def arrow_geometry(frm, to, bow, off0, off1, n=48) -> dict:
    """Polyline (trimmed at both node boxes) + arrowhead triangle."""
    p0 = (NODES[frm][1] + off0[0], NODES[frm][2] + off0[1])
    p2 = (NODES[to][1] + off1[0], NODES[to][2] + off1[1])
    dx, dy = p2[0] - p0[0], p2[1] - p0[1]
    d = math.hypot(dx, dy) or 1.0
    c = ((p0[0] + p2[0]) / 2 - bow * dy / d, (p0[1] + p2[1]) / 2 + bow * dx / d)
    pts = []
    for i in range(n + 1):
        t = i / n
        a, b, e = (1 - t) ** 2, 2 * (1 - t) * t, t ** 2
        pts.append((a * p0[0] + b * c[0] + e * p2[0],
                    a * p0[1] + b * c[1] + e * p2[1]))
    pts = [p for p in pts
           if not _inside(p, frm, NODE_PAD) and not _inside(p, to, NODE_PAD)]
    if len(pts) < 2:
        return {"line": [], "head": []}
    tip = pts[-1]
    # Direction at the tip, taken a few samples back for stability.
    ref = pts[max(0, len(pts) - 4)]
    ux, uy = tip[0] - ref[0], tip[1] - ref[1]
    u = math.hypot(ux, uy) or 1.0
    ux, uy = ux / u, uy / u
    base = (tip[0] - ux * HEAD_LEN, tip[1] - uy * HEAD_LEN)
    head = [tip,
            (base[0] - uy * HEAD_W, base[1] + ux * HEAD_W),
            (base[0] + uy * HEAD_W, base[1] - ux * HEAD_W)]
    # End the shaft at the arrowhead base so a thick line never pokes through.
    shaft = [p for p in pts
             if (p[0] - tip[0]) * ux + (p[1] - tip[1]) * uy <= -HEAD_LEN * 0.8]
    shaft.append(base)
    return {"line": shaft, "head": head}


def layout() -> dict:
    """Everything a renderer needs, JSON-serialisable."""
    r = lambda p: [round(p[0], 2), round(p[1], 2)]          # noqa: E731
    steps = []
    for st in STEPS:
        def geom(specs):
            out = []
            for a in specs:
                g = arrow_geometry(*a)
                out.append({"line": [r(p) for p in g["line"]],
                            "head": [r(p) for p in g["head"]]})
            return out
        # A reversible step also carries the same arrows drawn the other way.
        rev = [(to, frm, -bow, o1, o0) for frm, to, bow, o0, o1 in st["arrows"]] \
            if st.get("reversible") else []
        steps.append({"id": st["id"], "name": st["name"],
                      "pathway": st["pathway"], "arrows": geom(st["arrows"]),
                      "arrows_rev": geom(rev), "label": list(st["label"])})
    return {
        "xlim": list(XLIM), "ylim": list(YLIM), "line_step": LINE_STEP,
        "wide": (XLIM[1] - XLIM[0]) > 130,
        "nodes": [{"id": k, "label": v[0], "x": v[1], "y": v[2],
                   "w": v[3], "h": v[4]} for k, v in NODES.items()],
        "steps": steps,
    }
