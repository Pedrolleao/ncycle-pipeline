# common.smk — shared helpers used by both modes.

wildcard_constraints:
    sample=r"[A-Za-z0-9_.\-]+",


def sample_kind(sample):
    return SAMPLES[sample]["kind"]


def sample_path(sample):
    return SAMPLES[sample]["path"]


def sample_path2(sample):
    return SAMPLES[sample].get("path2", "")


def sample_is_protein_mode(sample):
    return sample_kind(sample) in ("protein", "nucleotide")


# Used by report.smk to dispatch the right calls-source per sample.
def calls_tsv_for(sample):
    return str(RESULTS / sample / "calls" / "ncycle_calls.tsv")
