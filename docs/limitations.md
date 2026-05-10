# Limitations

## Methodological

1. **MR p-value alone is insufficient.** Colocalization plus replication is required. AMBER / UNCLASSIFIED labels are a direct consequence of this rule.
2. **AF-mediated vs. off-axis is estimator-dependent.** `excess_off_axis_effect` depends on the accuracy of the AF → outcome MR signal; for an outcome where AF MR is weak (e.g. small HF subtype studies), the decomposition is unreliable.
3. **UKB sample overlap.** UKB-PPP discovery + UKB-containing outcome pairs require sensitivity analysis. A target that fails sensitivity is not Tier-1.
4. **Mixed-ancestry pQTL with EUR outcomes.** LD mismatch can break colocalization. The primary analysis is EUR-only; mixed-ancestry pQTL is used for sensitivity only.
5. **Trans-pQTL.** Pleiotropy risk is high; trans-pQTLs are excluded from the primary call.
6. **Nearest-gene assignment is not assumed.** Locus → gene mapping uses L2G + cis-QTL coloc + single-cell expression jointly.
7. **HFpEF GWAS power is limited.** HFpEF sample size is smaller than HFrEF; HFpEF harm signals are flagged "insufficient power" where appropriate.
8. **MEGASTROKE subtype power varies.** Small-vessel and large-artery subtypes are less well-powered than cardioembolic. Threshold-FDR can hide asymmetric statistical power across subtypes.

## Translational

1. **MR does not directly predict drug effect size.** "Lifelong genetic" exposure differs from "short-term drug" exposure in dose × duration.
2. **Direction may be correct, modality may not be.** Even when the analysis says "inhibit", small-molecule, antibody, ASO, and siRNA approaches can have very different tissue distributions.
3. **Off-target effects are not captured here.** Existing-drug side-effect profiles are kept in ChEMBL / DGIdb annotation but are not sufficient for de-novo drug design.
4. **Single-cell context alone is not causal evidence.** Causality is possible without expression, and expression itself does not exclude off-target risk.

## Data

1. **Liftover loss.** GRCh37 → GRCh38 liftover loses 1–2% of variants on average; > 5% triggers an audit.
2. **eQTLGen z-score → β conversion.** Some eQTLGen tables provide p and z but not β directly; harmonization must be careful.
3. **Bulk download is forbidden.** Disk-thrashing risk; streaming approaches (e.g. `phase2_proteome_wide_mr.py` style) are used at later phases.

## Reporting

1. **Class labels are threshold-dependent.** Changing thresholds changes labels; every report archives the threshold settings used.
2. **`therapeutic_balance_score` is for ranking, not classification.** Score comparisons are valid only when the same outcome panel was used.
