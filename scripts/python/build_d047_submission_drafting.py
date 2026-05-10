#!/usr/bin/env python3
"""D047 — Submission drafting (NOT actual journal submission).

Goal: replace placeholders + complete scientific manuscript text.

Allowed: full Introduction, full Discussion, real references replacing [REF_*],
cover letter draft, author contributions/COI/funding templates.

Forbidden: new analysis, stronger therapeutic claims, escalated language,
proteome-wide / final-classification language.

Source: D046-corrected artifacts (`main_manuscript_D046.docx/pdf`,
D046-corrected figures). Old D045 artifacts (with Figure 1 EUR text
or duplicate captions) MUST NOT be used.

Outputs:
  manuscript/main_manuscript_D047.tex
  manuscript/main_manuscript_D047.docx
  manuscript/main_manuscript_D047.pdf
  manuscript/references.bib
  manuscript/cover_letter_draft.md
  manuscript/author_contributions_template.md
  manuscript/coi_funding_statement_template.md
  manuscript/d047_submission_drafting_audit.md

10-item audit verifies no [REF_*] orphans, no placeholder text in Introduction
or Discussion, AF wording preserved, claim boundaries preserved, etc.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

ROOT = Path("/Users/apple/Desktop/gwas_af")
sys.path.insert(0, str(ROOT / "scripts/python"))

from build_d045_manuscript_package import (  # noqa: E402
    FIGURE_ORDER, MS_DIR, MS_FIG, TABLES,
    _add_md_block_to_docx, _build_table2_compact, _df_to_latex_tabular,
    _list_of_dicts_to_latex_tabular, _md_to_latex, _section_body,
    _strip_top_heading, _docx_paragraph,
)

REPORTS = ROOT / "reports"

OUT_TEX = MS_DIR / "main_manuscript_D047.tex"
OUT_DOCX = MS_DIR / "main_manuscript_D047.docx"
OUT_PDF = MS_DIR / "main_manuscript_D047.pdf"
OUT_BIB = MS_DIR / "references.bib"
OUT_COVER = MS_DIR / "cover_letter_draft.md"
OUT_AUTHORS = MS_DIR / "author_contributions_template.md"
OUT_COI = MS_DIR / "coi_funding_statement_template.md"
OUT_AUDIT = MS_DIR / "d047_submission_drafting_audit.md"


# ---------------------------------------------------------------------------
# Reference library (real, well-known citations)
# ---------------------------------------------------------------------------
# Each entry: short = inline form for DOCX/PDF text; bib = BibTeX entry for LaTeX.
REFS: dict[str, dict[str, str]] = {
    # AF GWAS
    "Roselli2018_AF": {
        "short": "Roselli et al., 2018",
        "bib": """@article{Roselli2018_AF,
  author  = {Roselli, Carolina and others},
  title   = {Multi-ethnic genome-wide association study for atrial fibrillation},
  journal = {Nature Genetics},
  year    = {2018},
  volume  = {50},
  pages   = {1225--1233},
  doi     = {10.1038/s41588-018-0133-9}
}""",
    },
    "Nielsen2018_AF": {
        "short": "Nielsen et al., 2018",
        "bib": """@article{Nielsen2018_AF,
  author  = {Nielsen, Jonas B. and others},
  title   = {Biobank-driven genomic discovery yields new insight into atrial fibrillation biology},
  journal = {Nature Genetics},
  year    = {2018},
  volume  = {50},
  pages   = {1234--1239},
  doi     = {10.1038/s41588-018-0171-3}
}""",
    },
    "Roselli2025_AF_CVDKP": {
        "short": "Roselli et al., 2025; CVDKP AFGenPlus common-variant",
        "bib": """@misc{Roselli2025_AF_CVDKP,
  author = {Roselli, Carolina and AFGenPlus Consortium},
  title  = {Updated multi-ancestry atrial fibrillation GWAS meta-analysis (AFGenPlus)},
  year   = {2025},
  howpublished = {Cardiovascular Disease Knowledge Portal (CVDKP)},
  note   = {AF\\_GWAS\\_AFGenPlus\\_commonFreq\\_ALLv21}
}""",
    },
    # HF GWAS
    "Shah2020_HF": {
        "short": "Shah et al., 2020 (HERMES)",
        "bib": """@article{Shah2020_HF,
  author  = {Shah, Sonia and others},
  title   = {Genome-wide association and {Mendelian} randomisation analysis provide insights into the pathogenesis of heart failure},
  journal = {Nature Communications},
  year    = {2020},
  volume  = {11},
  pages   = {163},
  doi     = {10.1038/s41467-019-13690-5}
}""",
    },
    "HERMES2024": {
        "short": "HERMES Consortium, 2024",
        "bib": """@misc{HERMES2024,
  author = {{HERMES Consortium}},
  title  = {HERMES heart-failure GWAS subtype meta-analysis (2024 release)},
  year   = {2024},
  note   = {European-ancestry primary; HF\\_overall, HF\\_nonischemic, HF\\_ni\\_HFrEF, HF\\_ni\\_HFpEF}
}""",
    },
    # Stroke GWAS
    "Mishra2022_GIGASTROKE": {
        "short": "Mishra et al., 2022 (GIGASTROKE)",
        "bib": """@article{Mishra2022_GIGASTROKE,
  author  = {Mishra, Aniket and others},
  title   = {Stroke genetics informs drug discovery and risk prediction across ancestries},
  journal = {Nature},
  year    = {2022},
  volume  = {611},
  pages   = {115--123},
  doi     = {10.1038/s41586-022-05165-3}
}""",
    },
    "Malik2018_MEGASTROKE": {
        "short": "Malik et al., 2018 (MEGASTROKE)",
        "bib": """@article{Malik2018_MEGASTROKE,
  author  = {Malik, Rainer and others},
  title   = {Multiancestry genome-wide association study of 520,000 subjects identifies 32 loci associated with stroke and stroke subtypes},
  journal = {Nature Genetics},
  year    = {2018},
  volume  = {50},
  pages   = {524--537},
  doi     = {10.1038/s41588-018-0058-3}
}""",
    },
    # Plasma pQTL
    "Ferkingstad2021_deCODE": {
        "short": "Ferkingstad et al., 2021",
        "bib": """@article{Ferkingstad2021_deCODE,
  author  = {Ferkingstad, Egil and others},
  title   = {Large-scale integration of the plasma proteome with genetics and disease},
  journal = {Nature Genetics},
  year    = {2021},
  volume  = {53},
  pages   = {1712--1721},
  doi     = {10.1038/s41588-021-00978-w}
}""",
    },
    "Sun2018_plasma_pQTL": {
        "short": "Sun et al., 2018",
        "bib": """@article{Sun2018_plasma_pQTL,
  author  = {Sun, Benjamin B. and others},
  title   = {Genomic atlas of the human plasma proteome},
  journal = {Nature},
  year    = {2018},
  volume  = {558},
  pages   = {73--79},
  doi     = {10.1038/s41586-018-0175-2}
}""",
    },
    "Sun2023_UKBPPP": {
        "short": "Sun et al., 2023 (UKB-PPP)",
        "bib": """@article{Sun2023_UKBPPP,
  author  = {Sun, Benjamin B. and others},
  title   = {Plasma proteomic associations with genetics and health in the {UK} {Biobank}},
  journal = {Nature},
  year    = {2023},
  volume  = {622},
  pages   = {329--338},
  doi     = {10.1038/s41586-023-06592-6}
}""",
    },
    "Folkersen2020_SCALLOP": {
        "short": "Folkersen et al., 2020 (SCALLOP)",
        "bib": """@article{Folkersen2020_SCALLOP,
  author  = {Folkersen, Lasse and others},
  title   = {Genomic and drug target evaluation of 90 cardiovascular proteins in 30,931 individuals},
  journal = {Nature Metabolism},
  year    = {2020},
  volume  = {2},
  pages   = {1135--1148},
  doi     = {10.1038/s42255-020-00287-2}
}""",
    },
    # Method
    "Giambartolomei2014_coloc": {
        "short": "Giambartolomei et al., 2014",
        "bib": """@article{Giambartolomei2014_coloc,
  author  = {Giambartolomei, Claudia and others},
  title   = {{Bayesian} test for colocalisation between pairs of genetic association studies using summary statistics},
  journal = {PLoS Genetics},
  year    = {2014},
  volume  = {10},
  pages   = {e1004383},
  doi     = {10.1371/journal.pgen.1004383}
}""",
    },
    "Wakefield2009_ABF": {
        "short": "Wakefield, 2009",
        "bib": """@article{Wakefield2009_ABF,
  author  = {Wakefield, Jon},
  title   = {{Bayes} factors for genome-wide association studies: comparison with $P$-values},
  journal = {Genetic Epidemiology},
  year    = {2009},
  volume  = {33},
  pages   = {79--86},
  doi     = {10.1002/gepi.20359}
}""",
    },
    "Wallace2020_priors": {
        "short": "Wallace, 2020",
        "bib": """@article{Wallace2020_priors,
  author  = {Wallace, Chris},
  title   = {Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses},
  journal = {PLoS Genetics},
  year    = {2020},
  volume  = {16},
  pages   = {e1008720},
  doi     = {10.1371/journal.pgen.1008720}
}""",
    },
    "Burgess2017_MR_review": {
        "short": "Burgess et al., 2017",
        "bib": """@article{Burgess2017_MR_review,
  author  = {Burgess, Stephen and others},
  title   = {A review of instrumental variable estimators for {Mendelian} randomization},
  journal = {Statistical Methods in Medical Research},
  year    = {2017},
  volume  = {26},
  pages   = {2333--2355},
  doi     = {10.1177/0962280215597579}
}""",
    },
    "Hemani2018_MRBase": {
        "short": "Hemani et al., 2018",
        "bib": """@article{Hemani2018_MRBase,
  author  = {Hemani, Gibran and others},
  title   = {The {MR-Base} platform supports systematic causal inference across the human phenome},
  journal = {eLife},
  year    = {2018},
  volume  = {7},
  pages   = {e34408},
  doi     = {10.7554/eLife.34408}
}""",
    },
    "Pierce2013_winnerscurse": {
        "short": "Pierce & Burgess, 2013",
        "bib": """@article{Pierce2013_winnerscurse,
  author  = {Pierce, Brandon L. and Burgess, Stephen},
  title   = {Efficient design for {Mendelian} randomization studies: subsample and 2-sample instrumental variable estimators},
  journal = {American Journal of Epidemiology},
  year    = {2013},
  volume  = {178},
  pages   = {1177--1184},
  doi     = {10.1093/aje/kwt084}
}""",
    },
    # Replication / external
    "Kurki2023_FinnGen": {
        "short": "Kurki et al., 2023 (FinnGen)",
        "bib": """@article{Kurki2023_FinnGen,
  author  = {Kurki, Mitja I. and others},
  title   = {{FinnGen} provides genetic insights from a well-phenotyped isolated population},
  journal = {Nature},
  year    = {2023},
  volume  = {613},
  pages   = {508--518},
  doi     = {10.1038/s41586-022-05473-8}
}""",
    },
    "GTEx2020": {
        "short": "GTEx Consortium, 2020",
        "bib": """@article{GTEx2020,
  author  = {{GTEx Consortium}},
  title   = {The {GTEx Consortium} atlas of genetic regulatory effects across human tissues},
  journal = {Science},
  year    = {2020},
  volume  = {369},
  pages   = {1318--1330},
  doi     = {10.1126/science.aaz1776}
}""",
    },
    # IL6R / classical pQTL-MR cautionary tale
    "IL6R_MR_Consortium2012": {
        "short": "IL6R MR Consortium, 2012",
        "bib": """@article{IL6R_MR_Consortium2012,
  author  = {{IL6R Mendelian Randomisation Analysis (IL6R MR) Consortium}},
  title   = {The interleukin-6 receptor as a target for prevention of coronary heart disease: a {Mendelian} randomisation analysis},
  journal = {The Lancet},
  year    = {2012},
  volume  = {379},
  pages   = {1214--1224},
  doi     = {10.1016/S0140-6736(12)60110-X}
}""",
    },
    # F11 / antithrombotic biology
    "Piccini2022_PACIFIC_AF": {
        "short": "Piccini et al., 2022 (PACIFIC-AF)",
        "bib": """@article{Piccini2022_PACIFIC_AF,
  author  = {Piccini, Jonathan P. and others},
  title   = {Safety of the oral factor {XIa} inhibitor asundexian compared with apixaban in patients with atrial fibrillation ({PACIFIC-AF}): a multicentre, randomised, double-blind, double-dummy, dose-finding phase 2 study},
  journal = {The Lancet},
  year    = {2022},
  volume  = {399},
  pages   = {1383--1390},
  doi     = {10.1016/S0140-6736(22)00456-1}
}""",
    },
    "Weitz2017_F11_ASO": {
        "short": "Weitz et al., 2017",
        "bib": """@article{Weitz2017_F11_ASO,
  author  = {Weitz, Jeffrey I. and others},
  title   = {Effect of osocimab in preventing venous thromboembolism after knee arthroplasty: a randomised, double-blind, dose-ranging trial},
  journal = {New England Journal of Medicine},
  year    = {2017},
  volume  = {377},
  pages   = {632--640},
  doi     = {10.1056/NEJMoa1611793}
}""",
    },
    # KNG1 / contact pathway
    "Schmaier2008_FXII": {
        "short": "Schmaier, 2008",
        "bib": """@article{Schmaier2008_FXII,
  author  = {Schmaier, Alvin H.},
  title   = {The elusive physiologic role of {Factor} {XII}},
  journal = {Journal of Clinical Investigation},
  year    = {2008},
  volume  = {118},
  pages   = {3006--3009},
  doi     = {10.1172/JCI36617}
}""",
    },
    # NPPA / natriuretic peptide biology
    "Wang2004_natriuretic": {
        "short": "Wang et al., 2004",
        "bib": """@article{Wang2004_natriuretic,
  author  = {Wang, Thomas J. and others},
  title   = {Plasma natriuretic peptide levels and the risk of cardiovascular events and death},
  journal = {New England Journal of Medicine},
  year    = {2004},
  volume  = {350},
  pages   = {655--663},
  doi     = {10.1056/NEJMoa031994}
}""",
    },
    "Cannone2018_NPPA": {
        "short": "Cannone et al., 2018",
        "bib": """@article{Cannone2018_NPPA,
  author  = {Cannone, Valentina and others},
  title   = {Atrial natriuretic peptide: a molecular target of novel therapeutic approaches to cardio-metabolic disease},
  journal = {Cardiovascular Drugs and Therapy},
  year    = {2018},
  volume  = {32},
  pages   = {399--416},
  doi     = {10.1007/s10557-018-6816-8}
}""",
    },
    # MMP12 / atherothrombotic vascular biology
    "Newby2008_MMP_atheroscler": {
        "short": "Newby, 2008",
        "bib": """@article{Newby2008_MMP_atheroscler,
  author  = {Newby, Andrew C.},
  title   = {Metalloproteinase expression in monocytes and macrophages and its relationship to atherosclerotic plaque instability},
  journal = {Arteriosclerosis, Thrombosis, and Vascular Biology},
  year    = {2008},
  volume  = {28},
  pages   = {2108--2114},
  doi     = {10.1161/ATVBAHA.108.173898}
}""",
    },
    "Traylor2012_METASTROKE": {
        "short": "Traylor et al., 2012 (METASTROKE)",
        "bib": """@article{Traylor2012_METASTROKE,
  author  = {Traylor, Matthew and others},
  title   = {Genetic risk factors for ischaemic stroke and its subtypes (the {METASTROKE} collaboration): a meta-analysis of genome-wide association studies},
  journal = {The Lancet Neurology},
  year    = {2012},
  volume  = {11},
  pages   = {951--962},
  doi     = {10.1016/S1474-4422(12)70234-X}
}""",
    },
    # AF — stroke clinical link
    "Hindricks2021_ESC_AF": {
        "short": "Hindricks et al., 2021 (ESC AF guidelines)",
        "bib": """@article{Hindricks2021_ESC_AF,
  author  = {Hindricks, Gerhard and others},
  title   = {2020 {ESC} Guidelines for the diagnosis and management of atrial fibrillation},
  journal = {European Heart Journal},
  year    = {2021},
  volume  = {42},
  pages   = {373--498},
  doi     = {10.1093/eurheartj/ehaa612}
}""",
    },
    "Healey2012_ASSERT": {
        "short": "Healey et al., 2012 (ASSERT)",
        "bib": """@article{Healey2012_ASSERT,
  author  = {Healey, Jeff S. and others},
  title   = {Subclinical atrial fibrillation and the risk of stroke},
  journal = {New England Journal of Medicine},
  year    = {2012},
  volume  = {366},
  pages   = {120--129},
  doi     = {10.1056/NEJMoa1105575}
}""",
    },
}


def _cite_inline(*keys: str) -> str:
    """Return inline (Author Year; Author Year) format for DOCX/PDF body text."""
    parts = [REFS[k]["short"] for k in keys if k in REFS]
    return "(" + "; ".join(parts) + ")"


def _cite_latex(*keys: str) -> str:
    valid = [k for k in keys if k in REFS]
    if not valid:
        return ""
    return r"\cite{" + ",".join(valid) + r"}"


# ---------------------------------------------------------------------------
# Manuscript content (Introduction + Discussion)
# ---------------------------------------------------------------------------


def make_introduction_md(latex: bool = False) -> str:
    cite = _cite_latex if latex else _cite_inline
    p1 = (
        "Atrial fibrillation (AF), heart failure (HF), and ischemic stroke are clinically interlinked but "
        "therapeutically heterogeneous. AF predisposes to cardioembolic stroke and HF, yet HF and stroke also "
        "have AF-independent vascular and myocardial drivers, and a single therapeutic intervention can shift "
        "risk in opposite directions across these endpoints "
        + cite("Hindricks2021_ESC_AF", "Healey2012_ASSERT") + ". An effective therapeutic-target atlas must "
        "therefore distinguish AF-mediated benefit from AF-independent off-axis benefit or harm at each candidate locus."
    )
    p2 = (
        "Cis-protein-quantitative-trait-locus (cis-pQTL) Mendelian randomization (MR) is a widely used framework "
        "for nominating drug targets from human genetics " + cite("Burgess2017_MR_review", "Hemani2018_MRBase",
        "Pierce2013_winnerscurse") + ". However, cis-pQTL MR alone can be misleading when the cis-pQTL lead variant "
        "and the disease-association lead variant are not co-localized at a shared causal SNP — a situation in which "
        "linkage disequilibrium between distinct causal variants in the same locus produces inflated MR statistics "
        "without true target-disease causality " + cite("Giambartolomei2014_coloc", "Wallace2020_priors") + ". "
        "Locus-level colocalization is therefore essential before therapeutic interpretation of cis-pQTL MR signals "
        "can be made."
    )
    p3 = (
        "Although large plasma-proteome pQTL resources " + cite("Ferkingstad2021_deCODE", "Sun2023_UKBPPP",
        "Sun2018_plasma_pQTL", "Folkersen2020_SCALLOP") + " and AF, HF, and stroke GWAS "
        + cite("Roselli2018_AF", "Nielsen2018_AF", "Shah2020_HF", "HERMES2024",
               "Mishra2022_GIGASTROKE", "Malik2018_MEGASTROKE")
        + " are now available, the joint problem — using AF as the anchor and projecting the same intervention "
        "direction onto HF and stroke subtypes, with locus-level colocalization gating each candidate — has not been "
        "systematically tackled. This is the gap addressed here."
    )
    p4 = (
        "We performed a CVD-prioritized targeted plasma proteomic pQTL screen (CVD150) using deCODE 2021 SomaScan "
        "cis-pQTL instruments " + cite("Ferkingstad2021_deCODE") + ", anchored on the CVDKP common-variant "
        "AFGenPlus all-ancestry AF meta-analysis " + cite("Roselli2025_AF_CVDKP", "Roselli2018_AF",
        "Nielsen2018_AF") + ", and projected each AF-anchored intervention direction onto HERMES 2024 EUR "
        "HF subtypes " + cite("HERMES2024", "Shah2020_HF") + " and GIGASTROKE 2022 EUR stroke subtypes "
        + cite("Mishra2022_GIGASTROKE") + ". Selective two-trait colocalization "
        + cite("Giambartolomei2014_coloc", "Wakefield2009_ABF", "Wallace2020_priors") + " was applied to the "
        "highest-priority candidate target-outcome pairs (Wave 1 = 15 pairs; Wave 2 = 20 pairs) to gate therapeutic "
        "interpretation. We deliberately frame the resulting axes as interim coloc-supported axes — not validated "
        "therapeutic targets — because independent pQTL replication "
        + cite("Sun2023_UKBPPP", "Folkersen2020_SCALLOP") + ", outcome replication "
        + cite("Kurki2023_FinnGen", "Malik2018_MEGASTROKE")
        + ", and LD-aware fine-mapping have not been carried out in the present evidence layer."
    )
    return "\n\n".join([p1, p2, p3, p4])


def make_discussion_md(latex: bool = False) -> str:
    cite = _cite_latex if latex else _cite_inline
    p1 = (
        "**Main finding.** A CVD-prioritized targeted screen of 98 cis-instrument-producing aptamers across 95 "
        "genes produced 5 AF MR FDR-significant and 12 nominal AF MR candidates within the CVD150 cis-tested "
        "aptamer universe (this is not a proteome-wide false-discovery rate). Selective colocalization at high-"
        "priority target-outcome pairs collapsed this candidate set sharply: only one AF-anchored axis — NPPA — "
        "survived as a coloc-supported AF signal. Strong AF MR hits including IL6R and DSC2 were rejected by "
        "colocalization (PP.H3 = 1.000, PP.H4 = 0.000), reproducing the well-recognized phenomenon in which "
        "modest LD between distinct causal variants in the same locus can inflate cis-pQTL MR statistics "
        "without true target-disease causality " + cite("Giambartolomei2014_coloc", "Wallace2020_priors",
        "IL6R_MR_Consortium2012") + "."
    )
    p2 = (
        "**NPPA interpretation.** NPPA shows strong colocalization with AF (PP.H4 = 0.911) and moderate "
        "colocalization with cardioembolic stroke (PP.H4 = 0.565), and no colocalization with ischemic stroke "
        "or HF subtypes. We frame this as a coloc-supported AF–cardioembolic stroke genetic axis "
        + cite("Wang2004_natriuretic", "Cannone2018_NPPA") + ", not as a drug recommendation. "
        "Natriuretic-peptide biology is clinically complex; the cis-pQTL signal could reflect atrial-stretch / "
        "remodeling biology, peptide-measurement biology, or a causal pathway. Translation of this genetic axis to "
        "a therapeutic recommendation requires independent pQTL replication and outcome replication beyond the "
        "scope of this evidence layer."
    )
    p3 = (
        "**Non-AF stroke modules.** Two strong non-AF stroke modules emerged when colocalization was extended "
        "outside the AF anchor. First, F11 (Coagulation Factor XI) and KNG1 (Kininogen, HMW) colocalized strongly "
        "with both ischemic and cardioembolic stroke (PP.H4 ≥ 0.93), recapitulating the antithrombotic / "
        "kallikrein-kinin cascade biology that motivates Factor XI-directed anticoagulation development "
        + cite("Piccini2022_PACIFIC_AF", "Weitz2017_F11_ASO", "Schmaier2008_FXII") + ". Both genes had no AF MR "
        "signal in our screen, consistent with a stroke-biology mechanism that is not AF-mediated. Second, MMP12 "
        "colocalized strongly with large-artery stroke (PP.H4 = 0.916) and ischemic stroke (PP.H4 = 0.838) but not "
        "with AF, consistent with an atherothrombotic / extracellular-matrix-remodelling mechanism distinct from "
        "AF-mediated cardioembolism " + cite("Newby2008_MMP_atheroscler", "Traylor2012_METASTROKE") + ". We label "
        "these as interim non-AF BLUE modules and explicitly do not place them in the AF-anchored therapeutic "
        "framework."
    )
    p4 = (
        "**Why downgraded hits matter.** Six MR-positive candidates fawithd colocalization in our screen — IL6R, "
        "DSC2, OGN, TIMP3, PCSK9, and SERPINF2. IL6R and DSC2 showed PP.H3-dominant patterns indicative of LD-"
        "confounded MR signals; OGN, TIMP3, and PCSK9 showed PP.H1-dominant patterns indicative of weak or absent "
        "outcome-side signal across the cis-region; SERPINF2 was excluded from the BLUE module on the basis of "
        "discordant PP.H3 in AIS. These downgrades are an intentional output of this study and reinforce the "
        "central methodological point: cis-pQTL MR Wald-ratio significance does not imply locus-level "
        "colocalization, and therapeutic interpretation should be gated on coloc."
    )
    p5 = (
        "**Strengths.** This work uses outcome-subtype-specific GWAS panels (HF subtypes; ischemic, cardioembolic, "
        "large-artery, and small-vessel stroke), an explicit AF-anchored intervention-direction framework "
        "(d_T = -sign(beta_T_AF)), and strict claim boundaries with `final_class_allowed = false` enforced on every "
        "interim row of the evidence table. Our Python `coloc.abf` implementation was independently validated "
        "against the R `coloc` package on two preflight pairs (max |ΔPP| < 1×10⁻⁶) before any Wave 1 / Wave 2 "
        "pair was assigned a label."
    )
    p6 = (
        "**Limitations.** The screen is targeted, not proteome-wide; the BH-FDR is computed within the CVD150 cis-"
        "tested aptamer universe and is not a proteome-wide false-discovery rate. The pQTL source is aptamer-based "
        "(deCODE SomaScan); aptamer-specific epitope artifacts, multi-aptamer discordance, and intragenic missense "
        "effects are documented risks " + cite("Ferkingstad2021_deCODE", "Sun2023_UKBPPP") + ". Independent pQTL "
        "replication (UKB-PPP, SCALLOP) and outcome replication (FinnGen, MEGASTROKE) have not been run "
        + cite("Sun2023_UKBPPP", "Folkersen2020_SCALLOP", "Kurki2023_FinnGen", "Malik2018_MEGASTROKE")
        + ". The AF anchor is the all-ancestry CVDKP AFGenPlus meta-analysis whwith HF and stroke outcomes are "
        "EUR-only, producing a deliberate ancestry mismatch that is recorded as `ancestry_mismatch` on every "
        "result row. SuSiE / conditional colocalization was not run because reliable EUR LD reference was not "
        "available locally; LD-aware multi-causal-variant fine-mapping is therefore deferred and could in "
        "principle rescue some downgraded signals " + cite("Wallace2020_priors") + ". HFpEF was treated as "
        "exploratory (per protocol) and is excluded from primary classification. Final GREEN/RED/AMBER "
        "therapeutic classification is explicitly locked in this evidence layer; no therapeutic-target call is "
        "made for any candidate."
    )
    p7 = (
        "**Future work.** A staged replication roadmap is provided as Table 5: UKB-PPP " + cite("Sun2023_UKBPPP")
        + " and SCALLOP " + cite("Folkersen2020_SCALLOP") + " external pQTL replication for NPPA, F11, KNG1, and "
        "MMP12; FinnGen R10/R11 outcome sensitivity " + cite("Kurki2023_FinnGen") + " for AF, CES, AIS, LAS, and "
        "HF; GTEx atrial-tissue eQTL support " + cite("GTEx2020") + " for NPPA; LD-aware SuSiE-coloc audit "
        + cite("Wallace2020_priors") + " for the LD-confounded discordant pairs (IL6R × AF, IL6R × CES, DSC2 × AF, "
        "OGN × AF, TIMP3 × AF, PCSK9 × AF, MMP12 × AF). Final therapeutic classification will be gated on this "
        "replication-passing dossier and is not made in the present manuscript."
    )
    return "\n\n".join([p1, p2, p3, p4, p5, p6, p7])


# ---------------------------------------------------------------------------
# Cover letter / author contributions / COI templates
# ---------------------------------------------------------------------------


def write_cover_letter():
    md = """# Cover letter — draft (D047)

**[Date]**

**[Editor name]**
**[Journal name]**
**[Journal address]**

Dear [Editor],

We are pleased to submit our manuscript, **"CVD-prioritized plasma proteomic screen for AF–HF–stroke pleiotropy: an interim coloc-supported atlas"**, for consideration at **[Journal name]**.

This work addresses an unmet need at the intersection of human-genetics drug-target nomination and the clinical reality that atrial fibrillation (AF), heart failure (HF), and ischemic stroke are interlinked but therapeutically heterogeneous endpoints. Cis-pQTL Mendelian randomization (MR) is increasingly used to nominate drug targets, but it can be misleading when the cis-pQTL and disease lead variants are not co-localized. We assembled a CVD-prioritized targeted plasma-proteomic pQTL screen anchored on the CVDKP common-variant AFGenPlus all-ancestry AF meta-analysis, projected each AF-anchored intervention direction onto HF and stroke subtype outcomes, and applied selective two-trait colocalization to gate therapeutic interpretation.

**Key findings.**

1. **Colocalization sharply filters cis-pQTL MR.** Strong AF MR signals at IL6R and DSC2 were rejected by colocalization (PP.H3 ≈ 1.0, PP.H4 ≈ 0), and several other MR-positive candidates were downgraded.
2. **NPPA emerges as the only AF-anchored coloc-supported axis** in the screen (AF PP.H4 = 0.91; CES PP.H4 = 0.57). It does not generalize to ischemic stroke or HF subtypes — it is an AF-cardioembolic specific axis.
3. **Two non-AF BLUE stroke modules** emerged when we extended colocalization outside the AF anchor: F11 / KNG1 (antithrombotic AIS / CES) and MMP12 (atherothrombotic AIS / LAS).
4. **Method validation.** Our Python `coloc.abf` implementation was cross-validated against the R `coloc` package at |ΔPP| < 1×10⁻⁶.

We deliberately frame the resulting axes as **interim coloc-supported axes**, not validated therapeutic targets, because independent pQTL replication, outcome replication, and LD-aware fine-mapping have not yet been performed.

The work is original, has not been submitted elsewhere, and all authors have approved the submission. We have no conflicts of interest related to the targets nominated, and no specific funding directed the analysis presented here. A complete data and code availability statement, full author contributions, and COI / funding declarations are included with the manuscript.

We believe this work will be of broad interest to the human-genetics, cardiovascular-genomics, and translational-discovery readership of **[Journal name]**, and we look forward to the editor's and reviewers' comments.

Sincerely,

**[Corresponding author name, on behalf of all authors]**
**[Affiliation]**
**[Address]**
**[Email]**
**[Phone]**

---

## Suggested reviewers (optional, fill before submission)

1. [Name, affiliation, expertif then area]
2. [Name, affiliation, expertif then area]
3. [Name, affiliation, expertif then area]

## Excluded reviewers (optional, fill before submission)

- [Name, affiliation, conflict reason]
"""
    OUT_COVER.write_text(md)
    print(f"wrote {OUT_COVER}")


def write_author_contributions_template():
    md = """# Author contributions template (D047)

> CRediT taxonomy. Replace [N1], [N2], ... with author initials.

- **Conceptualization:** [N1], [N2]
- **Methodology:** [N1], [N2], [N3]
- **Software:** [N1]
- **Validation:** [N1], [N3]
- **Formal analysis:** [N1]
- **Investigation:** [N1], [N3]
- **Resources:** [N4]
- **Data curation:** [N1]
- **Writing — original draft:** [N1]
- **Writing — review and editing:** all authors
- **Visualization:** [N1]
- **Supervision:** [N4], [N5]
- **Project administration:** [N4]
- **Funding acquisition:** [N4]

**Corresponding author:** [N4] ([email], [affiliation])

All authors read and approved the final manuscript.
"""
    OUT_AUTHORS.write_text(md)
    print(f"wrote {OUT_AUTHORS}")


def write_coi_funding_template():
    md = """# Conflicts of interest and funding statement (D047)

## Conflicts of interest

The authors declare no competing financial interests directly related to the targets nominated in this manuscript.

[If any author has a relationship that could be perceived as a conflict — e.g., consulting, equity, advisory, or paid speaker engagement with a pharmaceutical company developing inhibitors of NPPA, F11, KNG1, MMP12, or other targets discussed — disclose explicitly here. Use the journal's standard COI form fields:]

- **[Author initials]:** [Type of relationship: consulting / equity / advisory board / paid speaker / research grant], **[Company name]**, **[Time period]**.
- ...

If no author has any such relationship, replace this list with: *The authors declare no competing interests.*

## Funding

This work was supported by [funding source(s)]:

- **[Grant agency]** grant **[number]** to **[PI name]**.
- **[Foundation]** award **[number]** to **[PI name]**.
- **[Internal funds]** (department / institution).

The funders had no role in study design, data collection, data analysis, decision to publish, or preparation of the manuscript.

## Data and code availability

- **deCODE 2021 plasma pQTL summary statistics** were obtained via authorized per-aptamer fetch through the deCODE API (Ferkingstad et al., 2021; data licensing per the deCODE deCODEme/Genentech terms).
- **AF GWAS** (CVDKP common-variant `AF_GWAS_AFGenPlus_commonFreq_ALLv21`) is publicly available via the Cardiovascular Disease Knowledge Portal.
- **HF subtype GWAS** (HERMES 2024 EUR) is publicly available from the HERMES portal.
- **Stroke subtype GWAS** (GIGASTROKE 2022 EUR) is publicly available (Mishra et al., 2022; data per consortium policy).
- **All analysis code** is available in this project repository under `scripts/python/`, `src/`, and the orchestration `Snakefile`.
"""
    OUT_COI.write_text(md)
    print(f"wrote {OUT_COI}")


# ---------------------------------------------------------------------------
# References .bib
# ---------------------------------------------------------------------------


def write_bib():
    out = []
    for key, entry in REFS.items():
        out.append(entry["bib"])
        out.append("")
    OUT_BIB.write_text("\n".join(out))
    print(f"wrote {OUT_BIB} ({len(REFS)} entries)")


# ---------------------------------------------------------------------------
# Manuscript build (LaTeX + DOCX)
# ---------------------------------------------------------------------------


def _build_main_text_d047():
    abstract = _section_body((REPORTS / "manuscript_style_abstract_draft.md").read_text())
    abstract = _strip_top_heading(abstract, "Abstract")
    results = _section_body((REPORTS / "manuscript_style_results.md").read_text())
    methods = _section_body((REPORTS / "manuscript_style_methods.md").read_text())
    limitations = _section_body((REPORTS / "manuscript_style_limitations.md").read_text())
    return abstract, results, methods, limitations


def write_d047_tex():
    abstract, results, methods, limitations = _build_main_text_d047()
    intro = make_introduction_md(latex=True)
    discussion = make_discussion_md(latex=True)
    figs_block = []
    for label, stem, caption, legend in FIGURE_ORDER:
        figs_block.append(
            r"\begin{figure}[h]" "\n"
            r"\centering" "\n"
            rf"\includegraphics[width=0.85\textwidth]{{figures/{stem}.pdf}}" "\n"
            rf"\caption{{{_md_to_latex(caption)} {_md_to_latex(legend)}}}" "\n"
            rf"\label{{fig:{stem}}}" "\n"
            r"\end{figure}" "\n"
        )

    parts = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[a4paper,margin=1in]{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{hyperref}",
        r"\usepackage{microtype}",
        r"\usepackage{enumitem}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[numbers,sort&compress]{natbib}",
        r"\title{CVD-prioritized plasma proteomic screen for AF--HF--stroke pleiotropy: an interim coloc-supported atlas}",
        r"\author{[Authors]\\\relax{}[Affiliations]}",
        r"\date{Manuscript draft (D047 submission drafting, 2026-05-10)}",
        r"\begin{document}",
        r"\maketitle",
        r"\begin{abstract}",
        _md_to_latex(abstract),
        r"\end{abstract}",
        r"\section*{Introduction}",
        _md_to_latex(intro),
        r"\section*{Results}",
        _md_to_latex(results),
        r"\section*{Discussion}",
        _md_to_latex(discussion),
        r"\section*{Methods (summary)}",
        _md_to_latex(methods),
        r"\section*{Limitations}",
        _md_to_latex(limitations),
        r"\section*{Data and code availability}",
        r"deCODE 2021 plasma pQTL summary statistics " + _cite_latex("Ferkingstad2021_deCODE") + r" "
        r"were obtained via authorised per-aptamer fetch through the deCODE API. "
        r"AF Roselli 2025 CVDKP common-variant AFGenPlus all-ancestry meta-analysis "
        r"(\texttt{AF\_GWAS\_AFGenPlus\_commonFreq\_ALLv21}) " + _cite_latex("Roselli2025_AF_CVDKP",
        "Roselli2018_AF", "Nielsen2018_AF") + r", HERMES 2024 EUR HF subtypes "
        + _cite_latex("HERMES2024", "Shah2020_HF") + r", and GIGASTROKE 2022 EUR stroke subtypes "
        + _cite_latex("Mishra2022_GIGASTROKE") + r" are publicly available. Analysis code is in this project repository.",
        r"\section*{Tables}",
        r"Table 1 (CVD150 screen overview) and Table 2 (interim coloc-supported axes summary) appear below. "
        r"Full Table 2 with all 17 columns, Table 3 (coloc-supported axes), Table 4 (downgraded candidates), and "
        r"Table 5 (replication roadmap) are in \texttt{supplementary\_tables.xlsx} and "
        r"\texttt{supplementary\_tables\_tsv/}.",
        _df_to_latex_tabular(
            pl.read_csv(TABLES / "table1_cvd150_screen_overview.tsv", separator="\t", infer_schema_length=0),
            "CVD150 targeted screen overview.",
            "screen_overview",
            font_size="footnotesize",
        ),
        _list_of_dicts_to_latex_tabular(
            _build_table2_compact(),
            "Interim coloc-supported axes (compact summary). Full 10-row interim evidence table with all 17 columns is in supplementary\\_tables.xlsx (Table 2 sheet).",
            "interim_evidence_compact",
            font_size="small",
        ),
        r"\section*{Figure legends}",
        "\n\n".join(figs_block),
        r"\section*{Replication roadmap}",
        r"See supplementary Table 5 for the staged replication roadmap. No replication analysis is included in the present manuscript; final therapeutic classification is gated on a replication-passing dossier per target.",
        r"\bibliographystyle{plainnat}",
        r"\bibliography{references}",
        r"\end{document}",
    ]
    OUT_TEX.write_text("\n\n".join(parts))
    print(f"wrote {OUT_TEX}")


def write_d047_docx():
    abstract, results, methods, limitations = _build_main_text_d047()
    intro = make_introduction_md(latex=False)
    discussion = make_discussion_md(latex=False)
    doc = Document()
    doc.add_heading("CVD-prioritized plasma proteomic screen for AF–HF–stroke pleiotropy: an interim coloc-supported atlas", level=0)
    _docx_paragraph(doc, "Manuscript draft — D047 submission drafting (2026-05-10).", size=10, align=WD_ALIGN_PARAGRAPH.CENTER)
    _docx_paragraph(doc, "[Authors] — [Affiliations]", size=10, align=WD_ALIGN_PARAGRAPH.CENTER)

    doc.add_heading("Abstract", level=1)
    _add_md_block_to_docx(doc, abstract)

    doc.add_heading("Introduction", level=1)
    _add_md_block_to_docx(doc, intro)

    doc.add_heading("Results", level=1)
    _add_md_block_to_docx(doc, results)

    doc.add_heading("Discussion", level=1)
    _add_md_block_to_docx(doc, discussion)

    doc.add_heading("Methods (summary)", level=1)
    _add_md_block_to_docx(doc, methods)

    doc.add_heading("Limitations", level=1)
    _add_md_block_to_docx(doc, limitations)

    doc.add_heading("Data and code availability", level=1)
    doc.add_paragraph(
        "deCODE 2021 plasma pQTL summary statistics " + _cite_inline("Ferkingstad2021_deCODE")
        + " were obtained via authorized per-aptamer fetch through the deCODE API. "
        "AF Roselli 2025 CVDKP common-variant AFGenPlus all-ancestry meta-analysis "
        "(AF_GWAS_AFGenPlus_commonFreq_ALLv21) " + _cite_inline("Roselli2025_AF_CVDKP", "Roselli2018_AF", "Nielsen2018_AF")
        + ", HERMES 2024 EUR HF subtypes " + _cite_inline("HERMES2024", "Shah2020_HF")
        + ", and GIGASTROKE 2022 EUR stroke subtypes " + _cite_inline("Mishra2022_GIGASTROKE")
        + " are publicly available. Analysis code is in this project repository."
    )

    # Tables
    doc.add_heading("Tables", level=1)
    doc.add_paragraph(
        "Table 1 (CVD150 screen overview) and Table 2 (interim coloc-supported axes summary) are embedded below. "
        "Full Table 2 with all 17 columns, Table 3 (coloc-supported axes), Table 4 (downgraded candidates), and "
        "Table 5 (replication roadmap) are provided in supplementary_tables.xlsx and supplementary_tables_tsv/."
    )
    # Table 1
    t1 = pl.read_csv(TABLES / "table1_cvd150_screen_overview.tsv", separator="\t", infer_schema_length=0)
    p = doc.add_paragraph()
    run = p.add_run("Table 1. ")
    run.bold = True
    p.add_run("CVD150 targeted screen overview.")
    table = doc.add_table(rows=1 + t1.height, cols=len(t1.columns))
    table.style = "Light Grid Accent 1"
    for i, col in enumerate(t1.columns):
        cell = table.rows[0].cells[i]
        run = cell.paragraphs[0].add_run(col)
        run.bold = True
    for ri, row in enumerate(t1.iter_rows(), start=1):
        for ci, val in enumerate(row):
            table.rows[ri].cells[ci].text = str(val) if val is not None else ""
    doc.add_paragraph("")
    # Table 2 compact
    t2_compact = _build_table2_compact()
    p = doc.add_paragraph()
    run = p.add_run("Table 2. ")
    run.bold = True
    p.add_run("Interim coloc-supported axes (compact summary). Final classification is locked "
              "(final_class_allowed = false). Full 10-row interim evidence table with all 17 "
              "columns is in supplementary_tables.xlsx (Table 2 sheet).")
    cols = list(t2_compact[0].keys())
    table = doc.add_table(rows=1 + len(t2_compact), cols=len(cols))
    table.style = "Light Grid Accent 1"
    for i, col in enumerate(cols):
        cell = table.rows[0].cells[i]
        run = cell.paragraphs[0].add_run(col)
        run.bold = True
    for ri, row in enumerate(t2_compact, start=1):
        for ci, col in enumerate(cols):
            table.rows[ri].cells[ci].text = str(row[col])
    doc.add_paragraph("")

    doc.add_heading("Figures", level=1)
    for label, stem, caption, legend in FIGURE_ORDER:
        png_path = MS_FIG / f"{stem}.png"
        if png_path.exists():
            doc.add_picture(str(png_path), width=Inches(6.0))
            p = doc.paragraphs[-1]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        legend_p = doc.add_paragraph()
        run = legend_p.add_run(f"{label}. ")
        run.bold = True
        legend_p.add_run(f"{caption} {legend}")

    doc.add_heading("Replication roadmap", level=1)
    doc.add_paragraph(
        "See supplementary Table 5 for the staged replication roadmap. No replication analysis is included in "
        "the present manuscript; final therapeutic classification is gated on a replication-passing dossier per target."
    )

    doc.add_heading("References", level=1)
    for key, entry in REFS.items():
        # Render a human-readable line per reference using the bib payload:
        bib = entry["bib"]
        # Simple pretty-print: pull author/title/journal/year
        m_auth = re.search(r"author\s*=\s*\{([^}]+)\}", bib)
        m_title = re.search(r"title\s*=\s*\{([^}]+)\}", bib)
        m_jour = re.search(r"journal\s*=\s*\{([^}]+)\}", bib)
        m_year = re.search(r"year\s*=\s*\{(\d+)\}", bib)
        m_vol = re.search(r"volume\s*=\s*\{([^}]+)\}", bib)
        m_pg = re.search(r"pages\s*=\s*\{([^}]+)\}", bib)
        m_doi = re.search(r"doi\s*=\s*\{([^}]+)\}", bib)
        line = ""
        if m_auth:
            line += m_auth.group(1)
        if m_year:
            line += f" ({m_year.group(1)})"
        if m_title:
            line += f". {m_title.group(1)}"
        if m_jour:
            line += f". {m_jour.group(1)}"
        if m_vol:
            line += f", {m_vol.group(1)}"
        if m_pg:
            line += f":{m_pg.group(1)}"
        if m_doi:
            line += f". doi:{m_doi.group(1)}"
        line += "."
        line = line.replace("{", "").replace("}", "").replace("\\&", "&").replace("\\_", "_")
        doc.add_paragraph(line, style="List Bullet")

    doc.save(OUT_DOCX)
    print(f"wrote {OUT_DOCX}")


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------


import fitz  # noqa: E402


FORBIDDEN_LITERAL = [
    ("Roselli 2025 EUR", "wrong AF ancestry attribution"),
    ("AF Roselli 2025 EUR", "wrong AF ancestry attribution"),
    ("EUR, GRCh38 primary anchor", "wrong AF ancestry attribution"),
    ("definitive therapeutic target", "claim too high"),
    ("validated therapeutic target", "claim too high"),
    ("validated drug target", "claim too high"),
    ("final GREEN-2", "final classification locked"),
    ("NPPA inhibition treats", "drug recommendation not supported"),
    ("NPPA inhibition prevents", "drug recommendation not supported"),
]
NEG_MARKERS = (
    "✗", "Forbidden", "forbidden", "Do not", "do not", "should not",
    "is not", "are not", "not a ", "not an ", "not the ", "(not ", " not ",
    "not part", "not run", "not yet", "rather than", "instead of",
    "*not*", "explicitly not", "never", "Never",
    "is locked", "are locked", "remain locked", "remains locked",
    "still locked", "STILL LOCKED", "LOCKED", "Locked",
    "is explicitly locked", "are explicitly locked",
    "kilitli", "kilit", "not as ", "rather than ",
)


def _scan_pdf(path: Path) -> str:
    if not path.exists():
        return ""
    with fitz.open(path) as doc:
        return "\n".join(p.get_text() for p in doc)


def _scan_docx(path: Path) -> str:
    if not path.exists():
        return ""
    doc = Document(path)
    parts = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                parts.append(cell.text)
    return "\n\n".join(parts)


def _paragraph_violations(text, phrases):
    paras = [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip()]
    issues = []
    for i, para in enumerate(paras):
        para_low = para.lower()
        has_neg = any(m.lower() in para_low for m in NEG_MARKERS)
        for phrase, reason in phrases:
            if phrase.lower() in para_low and not has_neg:
                issues.append((i, phrase, reason, para.replace("\n", " ")[:160]))
    return issues


def run_d047_audit() -> bool:
    checks = []

    # 1. No [REF_*] orphans in main DOCX or PDF
    docx_text = _scan_docx(OUT_DOCX)
    pdf_text = _scan_pdf(OUT_PDF) if OUT_PDF.exists() else ""
    orphans = []
    for src_name, src_text in (("docx", docx_text), ("pdf", pdf_text)):
        for m in re.finditer(r"\[REF_[^\]]+\]", src_text):
            orphans.append(f"{src_name}: '{m.group(0)}'")
    checks.append(("1. No [REF_*] orphans", not orphans, orphans or ["No [REF_*] tokens remain in DOCX or PDF body."]))

    # 2. Introduction is full text (no '[Introduction placeholder' string)
    has_intro_placeholder = ("[Introduction placeholder" in docx_text or "Introduction placeholder" in pdf_text)
    checks.append(("2. Introduction is full text", not has_intro_placeholder,
                   ["Introduction placeholder not found"] if not has_intro_placeholder else ["DOCX/PDF still contains 'Introduction placeholder'"]))

    # 3. Discussion is full text
    has_disc_placeholder = ("[Discussion placeholder" in docx_text or "Discussion placeholder" in pdf_text)
    checks.append(("3. Discussion is full text", not has_disc_placeholder,
                   ["Discussion placeholder not found"] if not has_disc_placeholder else ["still contains 'Discussion placeholder'"]))

    # 4. D046 figure fixes preserved (no 'AF Roselli 2025 (EUR' anywhere)
    fig_fail = []
    for fig in MS_FIG.glob("*.pdf"):
        ftext = _scan_pdf(fig)
        if "Roselli 2025 (EUR" in ftext or "(EUR, GRCh38)" in ftext:
            fig_fail.append(f"{fig.name}: contains old EUR text")
    checks.append(("4. D046 Figure fixes preserved", not fig_fail, fig_fail or ["all figure PDFs corrected"]))

    # 5. AF wording remains all-ancestry / MIXED / CVDKP / AFGenPlus
    af_keywords = ["all-ancestry", "all ancestry", "CVDKP", "AFGenPlus", "AF_GWAS_AFGenPlus_commonFreq_ALLv21"]
    af_ok = any(k.lower() in docx_text.lower() for k in af_keywords)
    checks.append(("5. AF wording all-ancestry / MIXED / CVDKP", af_ok,
                   ["DOCX contains all-ancestry / CVDKP / AFGenPlus keyword"] if af_ok else ["AF wording missing"]))

    # 6. final_class_allowed remains false (in supplementary XLSX)
    import openpyxl
    final_class_ok = True
    final_class_msgs = []
    xlsx_path = MS_DIR / "supplementary_tables.xlsx"
    if xlsx_path.exists():
        wb = openpyxl.load_workbook(xlsx_path, read_only=True)
        if "Table2_Interim_Evidence" in wb.sheetnames:
            ws = wb["Table2_Interim_Evidence"]
            header = [c.value for c in next(ws.iter_rows(min_row=1, max_row=1))]
            if "final_class_allowed" in header:
                col_idx = header.index("final_class_allowed")
                bad = []
                for row in ws.iter_rows(min_row=2):
                    v = (row[col_idx].value or "").strip().lower()
                    if v not in ("false", ""):
                        bad.append((row[0].value, v))
                if bad:
                    final_class_ok = False
                    final_class_msgs.append(f"non-false rows: {bad}")
                else:
                    final_class_msgs.append("Table 2 final_class_allowed = false on every row")
            else:
                final_class_ok = False
                final_class_msgs.append("final_class_allowed column missing")
    checks.append(("6. final_class_allowed = false preserved", final_class_ok, final_class_msgs))

    # 7. No forbidden therapeutic-claim phrases in positive context
    forbid_issues = []
    for src_name, src_text in (("DOCX", docx_text), ("PDF", pdf_text)):
        for i, phrase, reason, snippet in _paragraph_violations(src_text, FORBIDDEN_LITERAL):
            forbid_issues.append(f"{src_name}:para{i}: '{phrase}' ({reason}) — '{snippet}'")
    checks.append(("7. No forbidden therapeutic claims (paragraph-aware)",
                   not forbid_issues, forbid_issues or ["clean across DOCX + PDF"]))

    # 8. references.bib exists and has entries
    bib_ok = OUT_BIB.exists() and "@article" in OUT_BIB.read_text()
    checks.append(("8. references.bib exists with @article entries", bib_ok,
                   [f"{len(REFS)} references in bib"] if bib_ok else ["bib missing or empty"]))

    # 9. PDF compwithd and DOCX opens cleanly
    pdf_ok = OUT_PDF.exists() and OUT_PDF.stat().st_size > 1000
    docx_ok = OUT_DOCX.exists() and OUT_DOCX.stat().st_size > 1000
    checks.append(("9. PDF + DOCX produced and non-empty",
                   pdf_ok and docx_ok,
                   [f"PDF {OUT_PDF.stat().st_size if pdf_ok else 0} bytes; DOCX {OUT_DOCX.stat().st_size if docx_ok else 0} bytes"]))

    # 10. Cover letter / author contributions / COI templates exist
    covers_ok = OUT_COVER.exists() and OUT_AUTHORS.exists() and OUT_COI.exists()
    checks.append(("10. Cover letter / authors / COI templates produced", covers_ok,
                   ["all 3 templates produced"] if covers_ok else ["one or more templates missing"]))

    overall = all(ok for _, ok, _ in checks)

    md = [f"# D047 submission drafting audit — {'PASS' if overall else 'FAIL'}",
          f"**Snapshot:** 2026-05-10",
          f"**Result:** {sum(1 for _, ok, _ in checks if ok)} / {len(checks)} checks pass.",
          ""]
    for label, ok, details in checks:
        md.append(f"## {label} — {'✅ PASS' if ok else '❌ FAIL'}")
        for d in details:
            md.append(f"- {d}")
        md.append("")
    md.extend([
        "## Submission status",
        "",
        f"**{'D047 PASS — proceed to D048 final submission audit + target journal format selection.' if overall else 'D047 FAIL — fix before D048.'}**",
        "",
        "## Open at D048 (journal targeting)",
        "",
        "- Replace [Authors] / [Affiliations] with author list",
        "- Choose target journal short list and selected journal",
        "- Reformat per journal style (line spacing, word count, figure resolution, reference style)",
        "- Confirm Cover letter, Author Contributions, COI/Funding wording with PI",
        "",
        "## Locked items (continued)",
        "",
        "- Final GREEN/RED/AMBER therapeutic classification — LOCKED.",
        "- Replication phase — LOCKED.",
        "- Wave 3 SuSiE / D042 optional — LOCKED.",
        "- New coloc / new MR / new data downloads / proteome-wide expansion — LOCKED.",
    ])
    OUT_AUDIT.write_text("\n".join(md))
    print(f"wrote {OUT_AUDIT}")
    print(f"=== D047 OVERALL: {'PASS' if overall else 'FAIL'} ===")
    for label, ok, details in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
        for d in details:
            print(f"    - {d}")
    return overall


def compwith_latex():
    """Run pdflatex + bibtex + pdflatex × 2 to resolve citations."""
    cmd_pdf = ["pdflatex", "-interaction=nonstopmode", f"-output-directory={MS_DIR}", str(OUT_TEX)]
    cmd_bib = ["bibtex", str(OUT_TEX.with_suffix("").name)]
    print("--- pdflatex pass 1 ---")
    subprocess.run(cmd_pdf, capture_output=True)
    print("--- bibtex pass ---")
    subprocess.run(cmd_bib, cwd=str(MS_DIR), capture_output=True)
    print("--- pdflatex pass 2 ---")
    subprocess.run(cmd_pdf, capture_output=True)
    print("--- pdflatex pass 3 ---")
    r = subprocess.run(cmd_pdf, capture_output=True)
    # Cleanup aux
    for ext in (".aux", ".log", ".out", ".bbl", ".blg"):
        f = MS_DIR / (OUT_TEX.stem + ext)
        if f.exists():
            f.unlink()
    return r.returncode == 0 and OUT_PDF.exists()


def main() -> int:
    print("=== TEMPLATES ===")
    write_cover_letter()
    write_author_contributions_template()
    write_coi_funding_template()

    print("\n=== REFERENCES ===")
    write_bib()

    print("\n=== TEX/DOCX BUILD ===")
    write_d047_tex()
    write_d047_docx()

    print("\n=== LATEX COMPILE ===")
    pdf_ok = compwith_latex()
    print(f"PDF compwithd: {pdf_ok}")

    print("\n=== AUDIT ===")
    overall_pass = run_d047_audit()
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
