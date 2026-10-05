"""Everything specific to the NITROGEN cycle that the figures and report need:
pathway order, palette, labels and output file names. The sulfur sister
pipeline has its own _domain.py and _cycle_model.py; every other figure /
report script is identical between the two.

Palette: the pathway hues are assigned in CAT_ORDER (the order the pathway
blocks sit next to each other) and were checked with the dataviz palette
validator in that order (adjacent CVD ΔE >= 8, normal-vision ΔE >= 15, both
modes). Three light-mode hues are below 3:1 on white, so colour never carries
identity alone: every coloured mark has a text label or a block header.
"""

CYCLE_LETTER = "N"
CYCLE_NAME = "nitrogen"
CALLS_TSV = "ncycle_calls.tsv"
LOCI_TSV = "ncycle_loci.tsv"

CAT_ORDER = [
    "nitrogen_fixation", "nitrification_ammonia", "nitrification_nitrite",
    "denitrification", "dnra", "anammox", "assimilatory",
    "ammonia_assimilation", "organic_n_mineralization",
]

PATHWAY_COLOR = {
    "nitrogen_fixation":        "#2a78d6",   # blue
    "nitrification_ammonia":    "#eb6834",   # orange
    "nitrification_nitrite":    "#1baf7a",   # aqua
    "denitrification":          "#eda100",   # yellow
    "dnra":                     "#e87ba4",   # magenta
    "anammox":                  "#008300",   # green
    "assimilatory":             "#4a3aa7",   # violet
    "ammonia_assimilation":     "#e34948",   # red
    "organic_n_mineralization": "#0093ad",   # cyan
}
# Same hues stepped for a dark surface (used by report.html only).
PATHWAY_COLOR_DARK = {
    "nitrogen_fixation":        "#3987e5",
    "nitrification_ammonia":    "#d95926",
    "nitrification_nitrite":    "#199e70",
    "denitrification":          "#c98500",
    "dnra":                     "#d55181",
    "anammox":                  "#008300",
    "assimilatory":             "#9085e9",
    "ammonia_assimilation":     "#e66767",
    "organic_n_mineralization": "#1aa3b5",
}
PATHWAY_LABEL = {
    "nitrogen_fixation":        "N fixation",
    "nitrification_ammonia":    "Ammonia oxidation",
    "nitrification_nitrite":    "Nitrite oxidation",
    "denitrification":          "Denitrification",
    "dnra":                     "DNRA",
    "anammox":                  "Anammox",
    "assimilatory":             "Assimilatory NO₃⁻/NO₂⁻ reduction",
    "ammonia_assimilation":     "Ammonia assimilation",
    "organic_n_mineralization": "Organic-N mineralization",
}

# Shorter form for the angled block headers of the overview grid.
PATHWAY_SHORT = dict(PATHWAY_LABEL, assimilatory="Assimilatory NO₃⁻/NO₂⁻")

# Tidy-ups applied to complex / module ids when they are shown as labels.
PRETTY_REPLACE = [("n2o", "N₂O"), ("dnra", "DNRA"), ("cnor", "cNOR"),
                  ("no reductase", "NO reductase")]

# Legend text for a module that is ruled out by an exclusion rule.
RULED_OUT_LABEL = "ruled out — an excluded gene is present"
