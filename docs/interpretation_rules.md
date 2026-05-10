# Target classification rules

For every target, `classify_target()` returns one of the labels below. The labels are **threshold-based** (FDR + PP.H4 + theta direction). The `therapeutic_balance_score` field is used only for within-class ranking, never for class assignment.

| Class | Definition |
|---|---|
| **GREEN-1_AF_BENEFIT_NEUTRAL_SAFETY** | Benefit on AF; no significant harm on HF or stroke; no signal flagged "unresolved". |
| **GREEN-2_SHARED_BENEFIT** | Benefit on AF + benefit on at least one HF or stroke subtype (FDR<0.05, PP.H4 ≥ 0.80). |
| **AMBER-1_AF_BENEFIT_UNRESOLVED** | Benefit on AF; HF or stroke significant at p<0.05 but coloc PP.H4 weak — undecided. |
| **RED-1_AF_BENEFIT_HF_HARM** | AF benefit + HF subtype harm (FDR<0.05, PP.H4 ≥ 0.50). |
| **RED-2_AF_BENEFIT_STROKE_HARM** | AF benefit + stroke subtype harm (FDR<0.05, PP.H4 ≥ 0.50). |
| **RED-3_AF_BENEFIT_MULTI_HARM** | AF benefit + harm on both HF and stroke subtypes. |
| **BLUE_NON_AF_SHARED_MEDIATOR** | No / weak AF signal but benefit on an HF or stroke subtype. Non-AF indication. |
| **PURPLE_HF_STROKE_BENEFIT_AF_HARM** | AF harm + HF or stroke benefit. Possible HF / stroke drug with an AF safety caveat. |
| **UNCLASSIFIED** | No signal strong enough; the reason is recorded in the `notes` field. |

## Clinical interpretation rules

### GREEN-2 → strongest candidates

- AF + HF + stroke shared benefit → "breakthrough candidate".
- If an existing drug or clinical-stage molecule exists, this becomes a repurposing opportunity.
- If the single-cell context is appropriate (atrial cardiomyocyte + fibroblast / endothelial), mechanistic support is strong.

### RED-1 → contractility / sarcomere / Ca²⁺ caution

- HFrEF harm suggests an off-axis effect via contractility or sarcomere biology.
- HFpEF harm suggests fibrosis or inflammation.
- High expression in ventricular cardiomyocytes (single-cell) strengthens the warning.

### RED-2 → vascular off-axis harm

- Cardioembolic stroke benefit is expected from AF reduction; small-vessel or large-artery harm signals an AF-independent vascular effect.
- Confirm with off-axis decomposition (`excess_off_axis_effect`).

### PURPLE → HF / stroke drug with AF safety flag

- Beneficial for HF or stroke, but AF risk increases → arrhythmia monitoring required during clinical development.

## Off-axis vs. AF-mediated interpretation

- Cardioembolic stroke + AF benefit: usually AF-mediated.
- Small-vessel / large-artery stroke + AF benefit: not AF-mediated; AF-independent vascular or coagulation effect.
- HFrEF harm + AF benefit: contractility channel.
- HFpEF harm + AF benefit: fibrosis / inflammation channel.

The master table carries an `excess_offaxis_<outcome>` column on every row; if it exceeds 0.05 (log OR), the row is flagged "off-axis dominant".
