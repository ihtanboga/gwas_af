# ARCHITECTURE — AFSHF-TARGETMAP

## Component map

```
┌──────────────┐    ┌──────────────────┐
│ config/*.yml │───▶│ pydantic schemas │ (src/afshf/config_schema.py)
└──────────────┘    └──────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 0: download.smk    │───▶│ data/registry/*.tsv     │
└──────────────────────────┘    └─────────────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 1A: harmonize.smk  │───▶│ data/processed/sumstats │ ← src/R/01_format_sumstats.R
│   src/afshf/harmonize.py │    │ (parquet, std schema)   │   uses src/afshf/harmonize.py rules
└──────────────────────────┘    └─────────────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 1B: instruments.smk│───▶│ data/processed/instr.   │ ← src/R/02_build_qtl_instruments.R
│   src/afshf/instruments  │    └─────────────────────────┘   (PLINK clump + F-stat)
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 2: overlap.smk     │───▶│ results/overlap/*.tsv   │ ← src/R/05_run_ldsc_or_hdl.R
│   LDSC / HDL / HyPrColoc │    └─────────────────────────┘   src/R/06_run_hyprcoloc.R
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 3: mr.smk          │───▶│ results/mr/             │ ← src/R/03_run_twosample_mr.R
│   TwoSampleMR            │    └─────────────────────────┘
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 4: coloc.smk       │───▶│ results/coloc/          │ ← src/R/04_run_coloc.R
│   coloc.abf + SuSiE      │    └─────────────────────────┘
└──────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────┐
│ Phase 5: annotation.smk                                  │
│   src/afshf/intervention_direction.py                    │
│   src/afshf/scoring.py  →  classify_target()             │
│   src/afshf/annotate_druggability.py                     │
│   src/afshf/annotate_singlecell.py                       │
│ Outputs:                                                 │
│   results/final/target_pleiotropy_master.tsv             │
│   results/final/red_flag_targets.tsv                     │
│   results/final/green_shared_benefit_targets.tsv         │
│   results/final/shared_loci_af_hf_stroke.tsv             │
│   results/final/druggability_annotation.tsv              │
└──────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────┐    ┌─────────────────────────┐
│ Phase 6: report.smk      │───▶│ reports/afshf_*.html    │
│   Quarto                 │    │ figures/*.pdf           │
└──────────────────────────┘    └─────────────────────────┘
```

## Standard summary-statistics schema

`src/afshf/harmonize.STANDARD_COLUMNS` and `src/R/utils.R::afshf_standard_columns` list the same 18 columns:

```
trait, source, variant_id, rsid, chr, pos, build,
effect_allele, other_allele, beta, se, or, pval, eaf,
n, n_cases, n_controls, ancestry
```

Every phase preserves this schema and may *add* columns (e.g. `f_stat`, `target_id`, `coloc_pph4`).

## Decision data flow (intervention direction)

```
beta_T_AF  ─────────► derive_action_direction() ───► d_T ∈ {-1, 0, +1}
                                                       │
                                                       ▼
              project_effect()  ◄──────────  d_T, beta_T_O
                       │
                       ▼
        TargetOutcomeEffect  → benefit_flag / harm_flag
                       │
        ┌──────────────┼──────────────┐
        ▼              ▼              ▼
   AF benefit?   shared benefit?    harm?      → aggregate_target_score()
        │              │              │              │
        └──────────────┼──────────────┘              ▼
                       ▼                    therapeutic_balance
                 classify_target()         (AF + Σshared - λΣharm)
                       │
                       ▼
        GREEN-1 / GREEN-2 / AMBER-1 / RED-1/2/3 / BLUE / PURPLE / UNCLASSIFIED
```

## Unit-test contract

The `tests/` folder commits to preserving these behaviours:

| Test | Assertion |
|---|---|
| Test 1 (spec §18) | beta_AF=+0.20 → d=-1, theta_AF=-0.20 (AF benefit), theta_HF=+0.10 (HF harm) |
| Test 2 (spec §18) | beta_AF=-0.30 → d=+1, shared benefit (AF + stroke) |
| Test 3 (spec §18) | A/T + EAF missing → drop; A/T + EAF=0.49 → drop; A/T + EAF=0.12 → keep |
| Test 4 (spec §18) | RED-1 classification (AF benefit + HF harm, FDR + PPH4 conditional) |
| Test 5 (spec §18) | GREEN-2 (AF + HF + stroke benefit, all FDR<0.05 and PPH4>0.8) |

Additional tests cover scoring components, palindromic edge cases, AMBER-1, RED-2/3, BLUE, PURPLE, UNCLASSIFIED.

## Checkpoints

- **`make test`** must produce all-green before pipeline launch.
- **After Phase 1**: `data/processed/sumstats/AF.parquet` has 18 columns and ≥10⁶ rows.
- **After Phase 3**: every (qtl_source, outcome) pair has `results/mr/*.tsv`; positive controls (PITX2, etc.) show significant AF signal.
- **After Phase 5**: `results/final/target_pleiotropy_master.tsv` has a populated `target_class` column for every target.
