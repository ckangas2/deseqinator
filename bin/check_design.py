#!/usr/bin/env python3
"""
Validate the samplesheet and the contrasts together, BEFORE any alignment,
and write the one design table every downstream step reads.

DESeq2 accepts almost any design and returns a results table. The failures
below all produce a plausible DE list rather than an error:

  1. Cell lines pooled into one `~ condition` model.
     Cell-line differences become residual variance (or, if unbalanced,
     leak into the condition effect). Each contrast is therefore scoped to
     one cell line, or to `*` with `cell_line` added to the design.

  2. Extra arms pulled into a contrast.
     GSE147507's A549-ACE2 series includes Ruxolitinib-treated samples.
     A contrast selects exactly two condition values; everything else in
     the samplesheet is excluded and the exclusion is written down.

  3. Too few replicates.
     n = 1 per arm runs (DESeq2 falls back on the trend) and gives
     confident-looking p-values. Below --min-replicates, stop.

  4. A covariate confounded with condition.
     batch == condition makes the condition term inestimable; DESeq2 then
     errors deep inside model fitting, or, with partial confounding,
     quietly attributes the effect to batch. Checked by matrix rank.

  5. Copy-paste errors: duplicate sample names, duplicate
     (cell_line, condition, replicate), mixed single/paired-end layouts.

Output: contrasts.validated.tsv (id, cell_line, reference, target, design,
covariates, samples, n_reference, n_target) and design_report.txt.
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
REQUIRED = ["sample", "fastq_1", "fastq_2", "cell_line", "condition", "replicate"]
CONTRAST_REQUIRED = ["id", "cell_line", "reference", "target"]
COVARIATES = ["batch"]          # optional samplesheet columns used as design terms


def read_csv(path):
    with open(path, newline="") as fh:
        rows = [{k.strip(): (v or "").strip() for k, v in r.items() if k is not None}
                for r in csv.DictReader(fh)]
        return rows


def design_rank_ok(rows, covariates, include_cell_line):
    """True if condition is estimable alongside the covariates.

    Builds the treatment-coded model matrix (intercept + covariate dummies +
    condition dummy) and checks that adding the condition column raises the
    rank. If it does not, condition lies in the span of the other terms and
    its effect cannot be separated from them.
    """
    cols = [np.ones(len(rows))]
    terms = list(covariates) + (["cell_line"] if include_cell_line else [])
    for term in terms:
        levels = sorted({r[term] for r in rows})
        for lv in levels[1:]:
            cols.append(np.array([1.0 if r[term] == lv else 0.0 for r in rows]))
    base = np.column_stack(cols)
    cond = np.array([1.0 if r["_arm"] == "target" else 0.0 for r in rows])
    full = np.column_stack([base, cond])
    return np.linalg.matrix_rank(full) > np.linalg.matrix_rank(base), full.shape[1], len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--contrasts", required=True)
    ap.add_argument("--min-replicates", type=int, default=3)
    ap.add_argument("--out-dir", default="design")
    args = ap.parse_args()

    problems = []
    notes = []

    # ------------------------------------------------------------ samplesheet
    ss = read_csv(args.samplesheet)
    if not ss:
        problems.append(f"samplesheet {args.samplesheet} has a header but no samples")
    else:
        missing_cols = [c for c in REQUIRED if c not in ss[0]]
        if missing_cols:
            problems.append(f"samplesheet is missing column(s): {missing_cols} "
                            f"(required: {','.join(REQUIRED)})")

    if problems:
        return report(problems, notes, args.out_dir)

    names = Counter(r["sample"] for r in ss)
    for n, c in names.items():
        if c > 1:
            problems.append(f"duplicate sample name: {n}")
        if not NAME_RE.match(n):
            problems.append(f"sample name '{n}' may only contain letters, digits, '.', '_' and '-'")
    for r in ss:
        for k in ("sample", "fastq_1", "cell_line", "condition", "replicate"):
            if not r[k]:
                problems.append(f"sample '{r['sample'] or '?'}' has an empty '{k}'")

    layouts = {bool(r["fastq_2"]) for r in ss}
    if len(layouts) > 1:
        pe = [r["sample"] for r in ss if r["fastq_2"]]
        se = [r["sample"] for r in ss if not r["fastq_2"]]
        problems.append(f"mixed single-end ({len(se)}) and paired-end ({len(pe)}) libraries. "
                        f"featureCounts counts fragments differently for each; split the run.")

    keys = Counter((r["cell_line"], r["condition"], r["replicate"]) for r in ss)
    for k, c in keys.items():
        if c > 1:
            problems.append(f"{c} samples share cell_line/condition/replicate {k}: "
                            f"probably a copy-paste error in the samplesheet")

    present_covariates = [c for c in COVARIATES if c in ss[0]]

    # ------------------------------------------------------------ contrasts
    cts = read_csv(args.contrasts)
    if not cts:
        problems.append(f"contrasts file {args.contrasts} has no contrasts")
        return report(problems, notes, args.out_dir)
    missing_cols = [c for c in CONTRAST_REQUIRED if c not in cts[0]]
    if missing_cols:
        problems.append(f"contrasts file is missing column(s): {missing_cols}")
        return report(problems, notes, args.out_dir)

    ids = Counter(c["id"] for c in cts)
    for i, n in ids.items():
        if n > 1:
            problems.append(f"duplicate contrast id: {i}")
        if not NAME_RE.match(i):
            problems.append(f"contrast id '{i}' may only contain letters, digits, '.', '_' and '-'")

    lines = sorted({r["cell_line"] for r in ss})
    validated = []
    used = set()
    for c in cts:
        cid, line, ref, tgt = c["id"], c["cell_line"], c["reference"], c["target"]
        if ref == tgt:
            problems.append(f"[{cid}] reference and target are both '{ref}'")
            continue
        if line != "*" and line not in lines:
            problems.append(f"[{cid}] cell_line '{line}' is not in the samplesheet "
                            f"(have: {lines})")
            continue
        pool = [r for r in ss if line == "*" or r["cell_line"] == line]
        conds = sorted({r["condition"] for r in pool})
        for lv in (ref, tgt):
            if lv not in conds:
                problems.append(f"[{cid}] condition '{lv}' not found for cell_line '{line}' "
                                f"(have: {conds}). Condition values are matched exactly.")
        sub = [dict(r, _arm="reference" if r["condition"] == ref else "target")
               for r in pool if r["condition"] in (ref, tgt)]
        n_ref = sum(r["_arm"] == "reference" for r in sub)
        n_tgt = len(sub) - n_ref
        if min(n_ref, n_tgt) < args.min_replicates:
            problems.append(f"[{cid}] {n_ref} reference vs {n_tgt} target samples; "
                            f"need >= {args.min_replicates} per arm (--min_replicates)")
            continue

        # Covariates enter the design only if they vary inside this contrast.
        covs = [v for v in present_covariates if len({r[v] for r in sub if r[v]}) > 1]
        for v in covs:
            if any(not r[v] for r in sub):
                problems.append(f"[{cid}] column '{v}' is blank for some samples in this contrast")
        multi_line = line == "*" and len({r["cell_line"] for r in sub}) > 1

        estimable, n_coef, n_samp = design_rank_ok(sub, covs, multi_line)
        terms = covs + (["cell_line"] if multi_line else []) + ["condition"]
        design = "~ " + " + ".join(terms)
        if not estimable:
            problems.append(f"[{cid}] condition is confounded with {terms[:-1]}: every "
                            f"difference between {ref} and {tgt} is also a difference in "
                            f"those terms, so the effect cannot be attributed to condition")
            continue
        if n_samp - n_coef < 2:
            problems.append(f"[{cid}] {design} uses {n_coef} coefficients on {n_samp} samples; "
                            f"too few residual degrees of freedom to estimate dispersion")
            continue

        excluded = [r["sample"] for r in ss if r not in pool or r["condition"] not in (ref, tgt)]
        if line != "*":
            dropped_here = sorted({r["condition"] for r in pool} - {ref, tgt})
            if dropped_here:
                notes.append(f"[{cid}] conditions in {line} not in this contrast "
                             f"(excluded): {dropped_here}")
        used.update(r["sample"] for r in sub)
        validated.append(dict(
            id=cid, cell_line=line, reference=ref, target=tgt, design=design,
            covariates=",".join(terms[:-1]),
            samples=",".join(r["sample"] for r in sub),
            reference_samples=",".join(r["sample"] for r in sub if r["_arm"] == "reference"),
            target_samples=",".join(r["sample"] for r in sub if r["_arm"] == "target"),
            n_reference=n_ref, n_target=n_tgt, n_excluded=len(excluded)))

    unused = [r["sample"] for r in ss if r["sample"] not in used]
    if unused:
        notes.append(f"{len(unused)} sample(s) are aligned and counted but enter no contrast: "
                     f"{unused}")

    return report(problems, notes, args.out_dir, validated)


def report(problems, notes, out_dir, validated=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lines = []
    for v in validated or []:
        lines.append(f"contrast {v['id']}: {v['target']} vs {v['reference']} "
                     f"in {v['cell_line']}   design {v['design']}")
        lines.append(f"  reference ({v['n_reference']}): {v['reference_samples']}")
        lines.append(f"  target    ({v['n_target']}): {v['target_samples']}")
    for n in notes:
        lines.append(f"note: {n}")
    for p in problems:
        lines.append(f"FAIL: {p}")
    lines.append(f"RESULT: {len(problems)} problem(s)")
    text = "\n".join(lines) + "\n"
    (out / "design_report.txt").write_text(text)
    sys.stdout.write(text)
    if problems:
        return 1
    cols = ["id", "cell_line", "reference", "target", "design", "covariates", "samples",
            "reference_samples", "target_samples", "n_reference", "n_target", "n_excluded"]
    with open(out / "contrasts.validated.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(validated)
    return 0


if __name__ == "__main__":
    sys.exit(main())
