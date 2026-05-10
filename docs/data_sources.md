# Data sources

| Trait / Tissue | Source | Access | Notes |
|---|---|---|---|
| AF | CVDKP — Roselli 2025 common-variant AF GWAS | https://kp4cd.org/ | Primary. UKB-inclusion risk → 2018 EUR fallback for sensitivity. |
| AF (sensitivity) | CVDKP — AF HRC 2018 EUR | kp4cd.org | n_case=55,114, n_ctrl=482,295. |
| AF (replication) | FinnGen R12 — endpoint I9_AF | https://finngen.gitbook.io/documentation/data-download | UKB-independent. |
| HF | CVDKP — HERMES 2024 trans-ancestry, with subtypes | kp4cd.org | n_case=153,174 (overall). HFrEF / HFpEF / non-ischemic subtypes in separate files. |
| HF (fallback) | CVDKP — HERMES 2020 EUR | kp4cd.org | n_case=47,309. |
| Stroke | CVDKP — MEGASTROKE 2018 EUR (AIS, CES, LAS, SVS) | kp4cd.org | Primary EUR. |
| Stroke (replication) | GIGASTROKE 2022 | kp4cd.org or original publication | MEGASTROKE follow-up. |
| pQTL — UKB-PPP | Synapse / UK Biobank RAP | `SYNAPSE_AUTH_TOKEN` env present | Discovery. n=54,219 EUR. UKB sample-overlap risk. |
| pQTL — deCODE | https://www.decode.com/summarydata/ | Form / authorization may be required | Replication, SomaScan v4. |
| pQTL — SCALLOP | http://www.scallop-consortium.com/ | Per cohort | Replication, Olink Inflammation / Cardio. |
| eQTL — eQTLGen | https://www.eqtlgen.org/ | Open | Whole blood, n=31,684. |
| eQTL — GTEx v8 | https://www.gtexportal.org/ | Open (full sumstats via dbGaP) | Heart_Atrial_Appendage is the most relevant tissue. |
| Drug annotation | Open Targets, ChEMBL, DGIdb, Pharos / TCRD | API + REST | Phase 5. |
| Single-cell | Heart Cell Atlas v2 (heartcellatlas.org), CELLxGENE | Open | 704k cells, conduction system included. |
| LD reference | 1000G phase 3 EUR PLINK bfile, optionally UKB EUR | Open (1000G), authorized (UKB) | Coloc + clumping. |

## License / citation conditions

### HERMES Consortium (HF, 2024 release)
- Access: open via the CVDKP `api.kpndataregistry.org` endpoint.
- **Publication condition:** "Journal publications using shared GWAS summary statistics must not precede publication of the main results manuscript." Any of our publications using these data must therefore appear after the HERMES manuscript.
- **Citation:** "HERMES consortium" banner authorship is required; named HERMES authors are at the discretion of the study team.

### MEGASTROKE 2018 (Malik et al., Nat Genet 2018, PMID 29531354) — primary stroke source
- Official URL: https://megastroke.org/download.html
- Access: terms-of-use checkbox approval (manual) — bypass not permitted.
- Expected files: `MEGASTROKE.1.AS.EUR.out`, `.2.AIS.EUR.out`, `.3.LAS.EUR.out`, `.4.CES.EUR.out`, `.5.SVS.EUR.out` (`.gz` accepted).

### ISGC distribution data-use language (verbatim, from `README_Traylor_2012.pdf`)

> "These data are provided 'as is', and without warranty, **for scientific and educational use only**. The use of these data for **commercial purposes is NOT allowed**. The inclusion of these data in **commercial databases is NOT allowed**. If you download these data, you acknowledge that these data will be used only for non-commercial research purposes; ... that **secondary distribution of the data without approval by secondary parties is prohibited**; and that the investigator will **cite the appropriate ISGC publication** in any communications or publications arising directly or indirectly from these data."
>
> Caution (sample overlap): "Some cases and controls were used by multiple studies across the ISGC. This is a **non-trivial, non-ignorable issue**."
>
> Citation: "(1) acknowledge that 'data were accessed through the **ISGC Cerebrovascular Disease Knowledge Portal**' and (2) cite the appropriate ISGC publication."

Pipeline application:
- All ISGC-derived data carry `sample_overlap_flag = possible_internal_ISGC_overlap`.
- Publications must cite "ISGC Cerebrovascular Disease Knowledge Portal" plus Malik 2018 (MEGASTROKE) or Traylor 2012 (METASTROKE).

### METASTROKE 2012 (Traylor et al., Lancet Neurology 2012) — not primary
- CVDKP fallback link: https://personal.broadinstitute.org/ryank/3490334.Traylor.2012.zip
- README: https://s3.amazonaws.com/broad-portal-resources/stroke/README_Traylor_2012.pdf — METASTROKE collaboration 2012 study.
- Status: not used as the primary stroke source; reviewed for data-use language only.

### Roselli 2025 / AFGen+ (AF)
- README "Publication: TBD" — citation will be added once the preprint / publication appears.
- CVDKP terms.

## Access tracking

`data/registry/access_log.md` records, for each dataset: download date, URL, n_total, build, license, and fwith checksum.

## Sample-overlap matrix

```
                AF_2025  AF_2018_EUR  HERMES_2024  HERMES_2020  MEGASTROKE  GIGASTROKE  FinnGen
UKB-PPP             ✓          ✗           ✓            ✓            ✗           ✓           ✗
deCODE              ✗          ✗           ✗            ✗            ✗           ✗           ✗
SCALLOP             minor      minor       minor        minor        minor       minor       ✗
eQTLGen             ✓          ✓           ✓            ✓            ✓           ✓           ✓ (some sub-cohorts)
GTEx                ✗          ✗           ✗            ✗            ✗           ✗           ✗
```

`✓` = serious sample overlap; sensitivity analysis is required. `✗` = no known overlap.
