# ncycle-pipeline make targets.
#
# `make regression` is the accuracy gate (scores the reference panel (test_panel)
# vs ground truth and FAILs below the floors in validation/test_regression.py).
# It is the real test for this pipeline.
#
# `make test_protein` is a lighter end-to-end SMOKE test: it runs the pipeline on
# the reference genomes downloaded by `make test_data` and confirms a non-empty
# results matrix — i.e. the pipeline runs start-to-finish. It does NOT check call
# accuracy; that is `make regression`.

TEST_DATA := ../test_data
REF_GENOMES := \
    $(TEST_DATA)/Aferrooxidans_ATCC23270.fna \
    $(TEST_DATA)/Mextorquens_AM1.fna \
    $(TEST_DATA)/Ecoli_K12.fna \
    $(TEST_DATA)/Soneidensis_MR1.fna \
    $(TEST_DATA)/Cmetallidurans_CH34.fna

NCBI_BASE := https://ftp.ncbi.nlm.nih.gov/genomes/all

URL_AFERR  := $(NCBI_BASE)/GCF/000/021/485/GCF_000021485.1_ASM2148v1/GCF_000021485.1_ASM2148v1_genomic.fna.gz
URL_MEXT   := $(NCBI_BASE)/GCF/000/022/685/GCF_000022685.1_ASM2268v1/GCF_000022685.1_ASM2268v1_genomic.fna.gz
URL_ECOLI  := $(NCBI_BASE)/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
URL_SHEW   := $(NCBI_BASE)/GCF/000/146/165/GCF_000146165.2_ASM14616v2/GCF_000146165.2_ASM14616v2_genomic.fna.gz
URL_CUPR   := $(NCBI_BASE)/GCF/000/196/015/GCF_000196015.1_ASM19601v1/GCF_000196015.1_ASM19601v1_genomic.fna.gz

PANEL := ../test_panel   # 33-genome reference panel (proteomes) for the regression gate

.PHONY: env test_data test_protein test regression regression-score kegg-contrast validate-panel verify-seeds clean clean_all

env:
	@mamba env create -f envs/ncycle.yaml 2>/dev/null || conda env create -f envs/ncycle.yaml

test_data: $(REF_GENOMES)

$(TEST_DATA)/Aferrooxidans_ATCC23270.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_AFERR) | gunzip > $@

$(TEST_DATA)/Mextorquens_AM1.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_MEXT) | gunzip > $@

$(TEST_DATA)/Ecoli_K12.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_ECOLI) | gunzip > $@

$(TEST_DATA)/Soneidensis_MR1.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_SHEW) | gunzip > $@

$(TEST_DATA)/Cmetallidurans_CH34.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_CUPR) | gunzip > $@

# End-to-end smoke test: run the pipeline on the reference genomes and confirm a
# non-empty results matrix (i.e. the pipeline runs start-to-finish). Accuracy is
# checked separately by `make regression`.
test_protein: test_data
	python run.py --input $(TEST_DATA) --mode protein --prodigal-mode single
	@python -c "import csv, sys, pathlib; \
m = pathlib.Path('results/multisample_matrix.tsv'); \
sys.exit('FAIL: results/multisample_matrix.tsv missing — pipeline did not finish') if not m.exists() else None; \
rows = list(csv.reader(m.open(), delimiter='\t')); \
sys.exit('FAIL: results matrix is empty') if len(rows) < 2 else None; \
print(f'OK: pipeline ran end-to-end — {len(rows)-1} data rows x {len(rows[0])-1} columns in results/multisample_matrix.tsv'); \
print('    (call accuracy is checked by: make regression)')"

# ── Seed-snapshot drift check (re-fetch UniProt + Pfam-A; diff against pin) ──
# Added 2026-05-30 (ROADMAP P5.1.5). Read-only; non-zero exit on drift. The
# pin lives at resources/.cache/seeds_manifest.tsv + pfam_release.txt and is
# bootstrapped via workflow/scripts/bootstrap_seed_cache.py. Builds prefer
# the cache (build_blast_db.py reads it before fetching), so silent UniProt
# revisions can no longer change the BLAST DB underneath us.
verify-seeds:
	python validation/verify_seeds.py

# ── Panel QC gate (per-genome identity check) ─────────────────────────────────
# Added 2026-05-30 (ROADMAP P5.0.6) after the 2026-05-30 audit caught a
# Burkholderia-thailandensis proteome under a Nwinogradskyi_Nb255.faa filename.
# Runs before any ground-truth regeneration so a mislabeled file fails fast.
validate-panel:
	python validation/validate_panel.py --panel $(PANEL)

# ── Accuracy regression gate (the real test for this pipeline) ────────────────
# Full gate: validate the panel, regenerate ground truth, run the pipeline on the
# reference panel, score, and FAIL (non-zero exit) if accuracy dropped below the
# floors in validation/test_regression.py. Needs the panel at $(PANEL) + the KOfam cache.
regression: validate-panel
	python validation/build_ground_truth.py
	python run.py --input $(PANEL) --skip-db-setup --cores 8
	python validation/score_ncycle.py
	python validation/test_regression.py

# Fast gate: re-score the EXISTING results/ matrix and check the floors (no
# pipeline run). Use after a scoring/ground-truth change when calls are current.
regression-score: validate-panel
	python validation/build_ground_truth.py
	python validation/score_ncycle.py
	python validation/test_regression.py

test: regression

# ── Dual-GT contrast: automated KEGG-derived GT vs the manual authoritative GT ──
# Builds the reproducible KEGG GT (32 of 39 panel genomes — the 7 anammox/comammox/
# Nitrospina/Nitrolancea/Methylacidiphilum/P.stutzeri-F2a specialists are absent from
# KEGG) and scores the pipeline against it. The manual ground_truth.tsv is the
# authoritative headline (make regression); this is the reproducible lower-bound
# contrast (KEGG over-calls homology traps, so trap precision deflates). Mirrors
# scycle-pipeline's automated-KEGG vs curated-function dual GT. Needs network.
kegg-contrast:
	python validation/build_kegg_ground_truth.py
	GT_FILE=kegg_ground_truth.tsv python validation/score_ncycle.py

clean:
	rm -rf results/*

clean_all: clean
	rm -rf resources/hmm/* resources/blast_db/* resources/nt_refs/*
