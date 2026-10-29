#!/usr/bin/env python3
"""
Pick the GSE147507 runs for the requested cell lines from ENA's run report,
and write runs.tsv plus a samplesheet.

Sample titles are parsed, never assumed. Every row is classified as
(series, cell line, condition, replicate) or rejected, and the selection is
refused unless each requested series holds exactly 3 Mock and 3 SARS-CoV-2
samples of one cell line. If GEO's titles differ from what this parser
expects, it stops and prints them rather than guessing.

    python3 scripts/select_runs.py ena_runs.tsv --series Series6,Series7 --list
"""
import argparse
import csv
import re
import sys

LINES = {
    "Calu3":     re.compile(r"calu[-_.]?3", re.I),
    "A549-ACE2": re.compile(r"a549[-_.]?ace2", re.I),
    "A549":      re.compile(r"a549(?![-_.]?ace2)", re.I),
    "NHBE":      re.compile(r"nhbe", re.I),
}
CONDS = {
    "Mock":       re.compile(r"\bmock\b|_mock_|mock[-_]", re.I),
    "SARS-CoV-2": re.compile(r"sars[-_]?cov[-_]?2", re.I),
}
EXCLUDE = re.compile(r"rux|ruxolitinib|iav|influenza|rsv|hpiv|ifn|ns1|ferret", re.I)
SERIES = re.compile(r"(series\d+)", re.I)
REP = re.compile(r"[_-](\d+)\s*$")


def classify(title):
    s = SERIES.search(title)
    series = s.group(1).capitalize() if s else ""
    if EXCLUDE.search(title):
        return series, None, None, None, "excluded treatment/virus/tissue"
    lines = [k for k, rx in LINES.items() if rx.search(title)]
    conds = [k for k, rx in CONDS.items() if rx.search(title)]
    rep = REP.search(title)
    if len(lines) != 1 or len(conds) != 1 or not rep:
        return series, None, None, None, f"unparsed (lines={lines}, conds={conds})"
    return series, lines[0], conds[0], rep.group(1), "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ena_tsv")
    ap.add_argument("--series", default="Series6,Series7",
                    help="GEO series labels to keep (Series7 = Calu3, Series6 = A549-ACE2)")
    ap.add_argument("--fastq-dir", default="fastq")
    ap.add_argument("--runs-out", default="config/runs.tsv")
    ap.add_argument("--samplesheet-out", default="samplesheet.csv")
    ap.add_argument("--downloads-out", default="config/downloads.tsv")
    ap.add_argument("--list", action="store_true", help="print every title and its parse, then stop")
    args = ap.parse_args()

    want = {s.strip().capitalize() for s in args.series.split(",") if s.strip()}
    with open(args.ena_tsv, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit("ENA report is empty")

    picked = []
    for r in rows:
        series, line, cond, rep, why = classify(r["sample_title"])
        if args.list:
            print(f"{r['run_accession']}\t{r['sample_title']}\t{series}\t{line}\t{cond}\t{rep}\t{why}")
        if series in want and why == "ok":
            picked.append(dict(r, series=series, cell_line=line, condition=cond, replicate=rep))
    if args.list:
        return 0

    problems = []
    for s in sorted(want):
        sub = [p for p in picked if p["series"] == s]
        lines = {p["cell_line"] for p in sub}
        n = {c: sum(p["condition"] == c for p in sub) for c in CONDS}
        if len(lines) != 1 or n != {"Mock": 3, "SARS-CoV-2": 3}:
            problems.append(f"{s}: expected one cell line with 3 Mock + 3 SARS-CoV-2, got "
                            f"lines={sorted(lines)} counts={n}")
    layouts = {p["library_layout"] for p in picked}
    if len(layouts) > 1:
        problems.append(f"mixed library layouts: {layouts}")
    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        print("Run with --list to see how every title was parsed.", file=sys.stderr)
        return 1

    picked.sort(key=lambda p: (p["cell_line"], p["condition"] != "Mock", p["replicate"]))
    with open(args.runs_out, "w") as fh:
        fh.write("sample\trun_accession\tsample_title\tseries\tcell_line\tcondition\treplicate\t"
                 "library_layout\tread_count\tfastq_bytes\n")
        for p in picked:
            p["sample"] = f"{p['cell_line']}_{'mock' if p['condition'] == 'Mock' else 'sars2'}_{p['replicate']}"
            fh.write("\t".join(str(p[k]) for k in (
                "sample", "run_accession", "sample_title", "series", "cell_line", "condition",
                "replicate", "library_layout", "read_count", "fastq_bytes")) + "\n")
    paired = layouts == {"PAIRED"}
    with open(args.downloads_out, "w") as fh:
        for p in picked:
            urls, md5s = p["fastq_ftp"].split(";"), p["fastq_md5"].split(";")
            if len(urls) != (2 if paired else 1):
                print(f"FAIL: {p['run_accession']} has {len(urls)} FASTQ file(s) for a "
                      f"{'paired' if paired else 'single'}-end layout", file=sys.stderr)
                return 1
            for u, m in zip(urls, md5s):
                fh.write(f"https://{u}\t{m}\t{args.fastq_dir}/{u.rsplit('/', 1)[1]}\n")
    with open(args.samplesheet_out, "w") as fh:
        fh.write("sample,fastq_1,fastq_2,cell_line,condition,replicate,batch,deposited_id\n")
        for p in picked:
            run = p["run_accession"]
            f1 = f"{args.fastq_dir}/{run}_1.fastq.gz" if paired else f"{args.fastq_dir}/{run}.fastq.gz"
            f2 = f"{args.fastq_dir}/{run}_2.fastq.gz" if paired else ""
            fh.write(f"{p['sample']},{f1},{f2},{p['cell_line']},{p['condition']},"
                     f"{p['replicate']},{p['series']},{p['sample_title']}\n")
    print(f"{len(picked)} runs ({'paired' if paired else 'single'}-end) -> "
          f"{args.runs_out}, {args.samplesheet_out}")
    for p in picked:
        print(f"  {p['sample']:<22} {p['run_accession']}  {p['sample_title']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
