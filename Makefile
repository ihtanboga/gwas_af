SHELL := /bin/bash
PY := .venv/bin/python
PIP := .venv/bin/pip
PYTEST := .venv/bin/pytest

.PHONY: help venv install test lint clean phase0 phase1 phase2 phase3 phase4 phase5 phase6 all

help:
	@echo "Common targets:"
	@echo "  venv          create local Python venv"
	@echo "  install       install Python deps inside the venv"
	@echo "  test          run pytest on src/afshf"
	@echo "  phase0        download metadata + dataset registry"
	@echo "  phase1        harmonize sumstats + build cis instruments"
	@echo "  phase2        outcome-outcome overlap (LDSC / HDL / HyPrColoc)"
	@echo "  phase3        target -> outcome MR"
	@echo "  phase4        QTL <-> outcome colocalization"
	@echo "  phase5        intervention direction + classification + scoring"
	@echo "  phase6        Quarto report"
	@echo "  all           run the entire pipeline"

venv:
	python3 -m venv .venv
	$(PIP) install --upgrade pip

install: venv
	$(PIP) install pytest pandas pydantic pyyaml numpy

test:
	$(PYTEST) tests/ -v

phase0:
	snakemake -j4 phase0_metadata

phase1:
	snakemake -j4 phase1_harmonize phase1_instruments

phase2:
	snakemake -j4 phase2_overlap

phase3:
	snakemake -j8 phase3_mr

phase4:
	snakemake -j8 phase4_coloc

phase5:
	snakemake -j4 phase5_classify

phase6:
	snakemake -j2 phase6_report

all:
	snakemake -j8

clean:
	rm -rf .pytest_cache __pycache__ src/afshf/__pycache__ tests/__pycache__
	find . -name "*.pyc" -delete
