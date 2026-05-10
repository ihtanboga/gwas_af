#!/usr/bin/env python3
"""D039 audit: re-run R coloc.abf properly (no sdY override) on the preflight
audit TSVs and compare against Python. Fixes JSON parsing by isolating the
last line of stdout that starts with '{'.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
AUDIT_DIR = ROOT / "results/coloc/preflight_audit"
PRIORS = {"p1": 1.0e-4, "p2": 1.0e-4, "p12": 1.0e-5, "W_pqtl": 0.0225, "W_outcome": 0.04}


def _logsumexp(xs):
    finite = [x for x in xs if x != float("-inf") and not math.isnan(x)]
    if not finite:
        return float("-inf")
    m = max(finite)
    return m + math.log(sum(math.exp(x - m) for x in finite))


def log_abf(z, V, W):
    if V <= 0:
        return float("-inf")
    r = W / (V + W)
    if r <= 0 or r >= 1:
        return float("-inf")
    return 0.5 * math.log(1.0 - r) + 0.5 * (z * z) * r


def python_coloc(rows, priors):
    n = len(rows)
    if n == 0:
        return {h: float("nan") for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")}
    log_abf_x, log_abf_y = [], []
    for r in rows:
        zx = r["beta"] / r["se"] if r["se"] > 0 else 0.0
        zy = r["beta_outcome"] / r["se_outcome"] if r["se_outcome"] > 0 else 0.0
        Vx = r["se"] ** 2; Vy = r["se_outcome"] ** 2
        log_abf_x.append(log_abf(zx, Vx, priors["W_pqtl"]))
        log_abf_y.append(log_abf(zy, Vy, priors["W_outcome"]))
    L_x = _logsumexp(log_abf_x)
    L_y = _logsumexp(log_abf_y)
    L_xy = _logsumexp([log_abf_x[i] + log_abf_y[i] for i in range(n)])
    log_sum = L_x + L_y
    if log_sum > L_xy and math.exp(L_xy - log_sum) < 1.0:
        L_distinct = log_sum + math.log(1.0 - math.exp(L_xy - log_sum))
    else:
        L_distinct = float("-inf")
    log_p1 = math.log(priors["p1"]); log_p2 = math.log(priors["p2"]); log_p12 = math.log(priors["p12"])
    log_terms = [0.0, log_p1 + L_x, log_p2 + L_y, log_p1 + log_p2 + L_distinct, log_p12 + L_xy]
    log_total = _logsumexp(log_terms)
    pp = [math.exp(x - log_total) for x in log_terms]
    return {"PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2], "PP.H3": pp[3], "PP.H4": pp[4]}


def r_coloc(rows, priors):
    rscript = shutil.which("Rscript")
    if not rscript:
        return None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        in_csv = td / "in.csv"
        with in_csv.open("w") as f:
            f.write("snp,beta_x,se_x,beta_y,se_y\n")
            for i, r in enumerate(rows):
                f.write(f"v{i},{r['beta']},{r['se']},{r['beta_outcome']},{r['se_outcome']}\n")
        # NO sdY: rely on R's default sd.prior=0.15 (quant) or 0.2 (cc).
        rcode = f"""
suppressMessages(library(coloc))
d <- read.csv("{in_csv}")
ds1 <- list(beta=d$beta_x, varbeta=d$se_x^2, type="quant", sdY=1, N=rep(35559, nrow(d)), snp=d$snp)
ds2 <- list(beta=d$beta_y, varbeta=d$se_y^2, type="cc", s=0.1, N=rep(100000, nrow(d)), snp=d$snp)
res <- suppressWarnings(suppressMessages(
  coloc.abf(ds1, ds2, p1={priors['p1']}, p2={priors['p2']}, p12={priors['p12']})
))
pp <- res$summary
out <- sprintf('JSON_OUT={{"PP.H0": %.8f, "PP.H1": %.8f, "PP.H2": %.8f, "PP.H3": %.8f, "PP.H4": %.8f}}',
               pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], pp["PP.H3.abf"], pp["PP.H4.abf"])
cat(out)
cat("\\n")
"""
        result = subprocess.run([rscript, "-e", rcode], capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"R FAIL: {result.stderr[:500]}", file=sys.stderr)
            return None
        for line in result.stdout.splitlines():
            if line.startswith("JSON_OUT="):
                return json.loads(line[len("JSON_OUT="):])
    return None


def main():
    files = sorted(AUDIT_DIR.glob("aligned_*.tsv"))
    print(f"Found {len(files)} audit files")
    out_rows = []
    for f in files:
        df = pl.read_csv(f, separator="\t")
        rows = df.to_dicts()
        py = python_coloc(rows, PRIORS)
        r = r_coloc(rows, PRIORS)
        diff = None
        passes = None
        if r is not None:
            diff = max(abs(py[h] - r[h]) for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))
            passes = diff < 1e-3
        tag = f.stem.replace("aligned_", "")
        out_rows.append({
            "pair": tag, "n_variants": len(rows),
            **{f"py_{h}": py[h] for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")},
            **({f"R_{h}": r[h] for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")} if r else
               {f"R_{h}": None for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")}),
            "max_abs_diff": diff, "passes_1e3": passes,
        })
        py_str = " ".join(f"{h}={py[h]:.4f}" for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))
        r_str = (" ".join(f"{h}={r[h]:.4f}" for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))
                 if r else "R unavailable")
        print(f"\n{tag} (n={len(rows)})")
        print(f"  Py: {py_str}")
        print(f"  R : {r_str}")
        if diff is not None:
            print(f"  max|Δ|={diff:.6f}  passes 1e-3? {passes}")

    out_path = ROOT / "results/coloc/decode_cvd150_wave1_preflight_validation.tsv"
    pl.DataFrame(out_rows).write_csv(out_path, separator="\t")
    print(f"\nwrote {out_path}")
    overall_pass = all(r["passes_1e3"] for r in out_rows if r["passes_1e3"] is not None)
    print(f"\nOverall: {'PASS' if overall_pass else 'FAIL'}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
