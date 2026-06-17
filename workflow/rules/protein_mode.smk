# protein_mode.smk — Mode A: protein input (.faa) or nucleotide input
# (.fna) → prodigal → .faa → hmmscan + DIAMOND-blastp fallback → apply_rules.

SCRIPTS = "workflow/scripts"


def proteome_for(sample):
    """Return the .faa path. For protein samples it is the input itself;
    for nucleotide samples it is the Prodigal output."""
    if sample_kind(sample) == "protein":
        return sample_path(sample)
    return str(RESULTS / sample / "prodigal" / f"{sample}.faa")


# ROADMAP P5.4.2: optionally query DIAMOND-blastp against only the HMM-hit
# proteins (config options.filter_blast_query_by_hmm) instead of the full
# proteome — ~50x fewer queries, no accuracy change. hmmscan + apply_rules
# (gene-coordinate synteny) always use the FULL proteome.
FILTER_BLAST_QUERY = bool(config.get("options", {}).get("filter_blast_query_by_hmm", False))


def blast_query_for(sample):
    """The FASTA used as the DIAMOND query: the HMM-filtered proteome when the
    toggle is on, else the full proteome."""
    if FILTER_BLAST_QUERY:
        return str(RESULTS / sample / "blast" / f"{sample}.hmm_filtered.faa")
    return proteome_for(sample)


rule prodigal:
    """Translate nucleotide input to predicted proteins."""
    input:
        lambda wc: sample_path(wc.sample),
    output:
        faa=RESULTS / "{sample}" / "prodigal" / "{sample}.faa",
        gff=RESULTS / "{sample}" / "prodigal" / "{sample}.gff",
    params:
        mode=lambda wc: SAMPLES[wc.sample].get("prodigal_mode", "single"),
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    log:
        RESULTS / "{sample}" / "prodigal" / "{sample}.log",
    shell:
        r"""
        prodigal -i {input} -a {output.faa} -f gff -o {output.gff} \
                 -p {params.mode} -q > {log} 2>&1
        """


rule hmmscan:
    """Scan a proteome against the concatenated N-cycle HMM database.
    Uses `hmmsearch` (profiles vs the proteome), not `hmmscan`: for a small
    profile DB scanned against thousands of proteins, hmmsearch is several-fold
    faster because it vectorises over the large sequence set per profile, while
    hmmscan pays profile-DB setup per query (the bottleneck at GTDB/MAG scale).
    The bitscores are identical, so calls are unchanged; `parse_hmmscan_domtbl`
    reads the swapped hmmsearch domtbl columns. hmmsearch reads the plain .hmm
    (no hmmpress needed). Output: domain-table (`--domtblout`)."""
    input:
        faa=lambda wc: proteome_for(wc.sample),
        hmm=config["paths"]["hmm_db"],
    output:
        tbl=RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.tsv",
    params:
        evalue=config["thresholds"]["hmmer"]["evalue"],
        dom_e=config["thresholds"]["hmmer"]["dom_evalue"],
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    threads: 4
    log:
        RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.log",
    shell:
        r"""
        hmmsearch --cpu {threads} -E {params.evalue} --domE {params.dom_e} \
                  --domtblout {output.tbl} {input.hmm} {input.faa} > {log} 2>&1
        """


rule filter_proteome_by_hmm:
    """ROADMAP P5.4.2 — write a proteome FASTA containing only proteins with
    >=1 hmmscan hit, used as the DIAMOND query (when options.filter_blast_query_by_hmm).
    Safe: every BLAST-db target has an HMM profile, so its true proteins are kept
    (the dropped >95% have zero HMM hits → zero target-specific BLAST evidence)."""
    input:
        faa=lambda wc: proteome_for(wc.sample),
        hmm=RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.tsv",
    output:
        faa=RESULTS / "{sample}" / "blast" / "{sample}.hmm_filtered.faa",
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    log:
        RESULTS / "{sample}" / "blast" / "{sample}.hmm_filtered.log",
    shell:
        r"""
        python {SCRIPTS}/filter_proteome_by_hmm.py \
            --proteome {input.faa} --hmmscan {input.hmm} \
            --out {output.faa} > {log} 2>&1
        """


rule diamond_blastp_unstable:
    """DIAMOND-blastp the proteome against the unstable_refs database
    (9 weak/⚠️ Pfam targets)."""
    input:
        faa=lambda wc: blast_query_for(wc.sample),
        db=Path(config["paths"]["blast_unstable"]).with_suffix(".dmnd"),
    output:
        tsv=RESULTS / "{sample}" / "blast" / "{sample}.unstable.tsv",
    params:
        evalue=config["thresholds"]["blast"]["evalue"],
        qcov=config["thresholds"]["blast"]["qcov_hsp_perc"],
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    threads: 4
    log:
        RESULTS / "{sample}" / "blast" / "{sample}.unstable.log",
    shell:
        r"""
        diamond blastp --quiet --threads {threads} --evalue {params.evalue} \
                       --query-cover {params.qcov} \
                       --db {input.db} --query {input.faa} \
                       --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                                   qstart qend sstart send evalue bitscore \
                       --out {output.tsv} 2> {log}
        """


rule diamond_blastp_gated:
    """DIAMOND-blastp the proteome against the blast_gated_refs database
    (no-Pfam targets + targets with `requires_blast_for_confirmation: true`).
    Permissive pident threshold applied downstream in apply_rules.py."""
    input:
        faa=lambda wc: blast_query_for(wc.sample),
        db=Path(config["paths"]["blast_gated"]).with_suffix(".dmnd"),
    output:
        tsv=RESULTS / "{sample}" / "blast" / "{sample}.blast_gated.tsv",
    params:
        evalue=config["thresholds"]["blast"]["evalue"],
        qcov=config["thresholds"]["blast"]["qcov_hsp_perc"],
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    threads: 4
    log:
        RESULTS / "{sample}" / "blast" / "{sample}.blast_gated.log",
    shell:
        r"""
        diamond blastp --quiet --threads {threads} --evalue {params.evalue} \
                       --query-cover {params.qcov} \
                       --db {input.db} --query {input.faa} \
                       --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                                   qstart qend sstart send evalue bitscore \
                       --out {output.tsv} 2> {log}
        """


rule apply_rules_protein:
    """Combine hmmscan + DIAMOND-blastp evidence and assign per-target status."""
    input:
        hmm=RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.tsv",
        blast_u=RESULTS / "{sample}" / "blast" / "{sample}.unstable.tsv",
        blast_g=RESULTS / "{sample}" / "blast" / "{sample}.blast_gated.tsv",
        targets=config["paths"]["targets_yaml"],
        # The proteome .faa doubles as the gene-coordinate source: for nucleotide
        # samples it is the Prodigal output (coords in the headers → nxr/nar synteny
        # resolution); for protein samples it is the pre-called input (no coords).
        faa=lambda wc: proteome_for(wc.sample),
    output:
        calls=RESULTS / "{sample}" / "calls" / "ncycle_calls.tsv",
    params:
        sample="{sample}",
        mode="protein",
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    log:
        RESULTS / "{sample}" / "calls" / "apply_rules.log",
    shell:
        r"""
        python {SCRIPTS}/apply_rules.py \
            --sample {params.sample} --mode {params.mode} \
            --hmm {input.hmm} \
            --blast-unstable {input.blast_u} \
            --blast-gated {input.blast_g} \
            --targets {input.targets} \
            --gene-coords {input.faa} \
            --out {output.calls} > {log} 2>&1
        """


