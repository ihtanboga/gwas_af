# Methods

## Standard summary-statistics schema

`src/afshf/harmonize.STANDARD_COLUMNS` lists 18 columns. All harmonized files are stored as parquet on GRCh38.

## Intervention direction

```
beta_T_AF = MR beta for genetically increased target -> AF log(OR)
d_T       = -sign(beta_T_AF)
```

- `d_T = -1`: drug intervention should inhibit / decrease the target
- `d_T = +1`: drug intervention should activate / increase the target
- `d_T = 0`: AF signal too weak (p ≥ 0.05) — target is not currently actionable

For each outcome O:

```
theta_T_O = d_T * beta_T_O
OR_T_O    = exp(theta_T_O)
```

Interpretation:
- `theta_T_AF < 0` → benefit on AF
- `theta_T_HF > 0` → harm on HF
- `theta_T_stroke > 0` → harm on stroke
- `theta_T_O < 0` → benefit on outcome O

## Off-axis decomposition

```
expected_AF_mediated  = theta_T_AF * beta_AF_to_outcome
excess_off_axis       = theta_T_outcome - expected_AF_mediated
```

- `excess > 0` → outcome harm exceeds the AF-mediated expectation → off-axis harm.
- `excess ≈ 0` → outcome effect is AF-mediated; the target appears safe.
- `excess < 0` → outcome benefit exceeds the AF-mediated expectation → strong shared-benefit candidate.

## Colocalization gating

- Strong: PP.H4 ≥ 0.80
- Moderate: PP.H4 ≥ 0.50
- LD-confounding flag: PP.H3 ≥ 0.50 and PP.H4 < 0.50

GREEN-2 classification requires PP.H4 ≥ 0.80 on shared-benefit outcomes. RED classification accepts PP.H4 ≥ 0.50 on harm outcomes (the lower bar avoids missing harm signals).

## Mendelian randomization

- Primary: single-SNP → Wald ratio; multi-SNP → IVW random-effects.
- Sensitivity: weighted median, MR-Egger, weighted mode, leave-one-out, heterogeneity Q, Steiger (flag-only), MR-PRESSO (≥4 SNPs).
- BH-FDR per outcome.
- F-statistic > 10 required.

## Colocalization

- Primary: `coloc.abf` with `p1 = p2 = 1×10⁻⁴`, `p12 = 5×10⁻⁶`.
- Sensitivity: `p12 ∈ {1×10⁻⁶, 5×10⁻⁶, 1×10⁻⁵, 1×10⁻⁴}`.
- If LD reference is available, `coloc.susie` is run as a secondary analysis.
- All SNPs in the regional window are included; significance-based filtering inside the coloc region is forbidden.

## Genetic correlation and multi-trait sharing

- LDSC or HDL → outcome × outcome rg matrix.
- HyPrColoc → multi-trait locus sharing.
- Per-locus lead-variant overlap and direction-concordance table.
