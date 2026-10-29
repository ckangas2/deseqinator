#!/usr/bin/env python3
"""
Assert the outcome of `-profile test` against the planted truth.

    python3 tests/check_test_results.py results

Every assertion is about what the pipeline found, not that it ran: planted
DE genes recovered with the right sign, few false positives, the planted
gene sets enriched in the planted direction, the unchanged set not, decoy
samples kept out of the contrast, viral load where it was put.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRAST = "lineA_inf_vs_mock"


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(res):
    res = Path(res)
    fails = []
    def check(cond, msg):
        print(("  PASS  " if cond else "  FAIL  ") + msg)
        if not cond:
            fails.append(msg)

    print("design")
    v = tsv(res / "design/contrasts.validated.tsv")
    want = {f"lineA_{c}_{r}" for c in ("mock", "inf") for r in (1, 2, 3)}
    check(len(v) == 1 and set(v[0]["samples"].split(",")) == want,
          "contrast holds exactly the six lineA mock/infected samples (no Rux, no lineB)")
    check(v[0]["design"] == "~ condition", f"design is '~ condition' (got '{v[0]['design']}')")

    print("strandedness")
    s = tsv(res / "qc/strandedness.tsv")
    check(len(s) == 9 and all(r["call"] == "reverse" for r in s), "all 9 libraries called reverse")

    print("viral load")
    vl = {r["sample"]: float(r["viral_frac"]) for r in tsv(res / "qc/viral_load.tsv")}
    mock = [k for k in vl if "mock" in k]
    inf = [k for k in vl if "mock" not in k]
    check(all(vl[k] == 0 for k in mock), "mock libraries carry no viral reads")
    check(all(0.15 < vl[k] < 0.25 for k in inf), "infected libraries are ~20% viral")

    print("counting")
    a = tsv(res / "counts/assignment.tsv")
    check(all(float(r["assigned_frac"]) > 0.6 for r in a), "assignment > 60% in every library")

    print("differential expression")
    truth = {r["gene_id"]: r["planted"] for r in tsv(ROOT / "tests/data/truth.tsv")}
    de = tsv(res / f"deseq2/{CONTRAST}/{CONTRAST}.deseq2.tsv")
    sig = {r["gene_id"]: float(r["log2FoldChange"]) for r in de
           if r["padj"] != "NA" and float(r["padj"]) < 0.05}
    hit_up = sum(1 for g, p in truth.items() if p == "up" and sig.get(g, 0) > 1)
    hit_dn = sum(1 for g, p in truth.items() if p == "down" and sig.get(g, 0) < -1)
    fp = sum(1 for g in sig if truth.get(g) == "none" and abs(sig[g]) >= 1)
    check(hit_up >= 17, f"{hit_up}/20 planted up genes recovered (padj<0.05, LFC>1)")
    check(hit_dn >= 17, f"{hit_dn}/20 planted down genes recovered (padj<0.05, LFC<-1)")
    check(fp <= 4, f"{fp} false positives among 160 unchanged genes (|LFC|>=1)")
    sm = tsv(res / f"deseq2/{CONTRAST}/{CONTRAST}.deseq2_summary.tsv")[0]
    check(int(sm["label_swaps"]) == 9, f"9 label swaps evaluated (got {sm['label_swaps']})")
    check(int(sm["n_up"]) + int(sm["n_down"]) > float(sm["swap_max_de"]),
          f"real DE ({int(sm['n_up']) + int(sm['n_down'])}) exceeds every label swap "
          f"(max {sm['swap_max_de']})")

    print("gene sets")
    fg = {r["pathway"]: r for r in tsv(res / f"fgsea/{CONTRAST}/{CONTRAST}.fgsea.tsv")}
    def p(k): return float(fg[k]["padj"])
    def nes(k): return float(fg[k]["NES"])
    check(p("HALLMARK_TEST_UP") < 0.05 and nes("HALLMARK_TEST_UP") > 0, "UP set enriched, NES > 0")
    check(p("HALLMARK_TEST_DOWN") < 0.05 and nes("HALLMARK_TEST_DOWN") < 0, "DOWN set enriched, NES < 0")
    check(p("HALLMARK_TEST_NULL") >= 0.05, "unchanged set not enriched")
    fs = tsv(res / f"fgsea/{CONTRAST}/{CONTRAST}.fgsea_summary.tsv")[0]
    check(int(fs["rank_shuffle_hits"]) == 0, "rank-shuffle control: no significant sets")

    print("reconciliation")
    rc = tsv(res / "qc/reconciliation.tsv")
    check(all(r["verdict"] in ("ok", "low") for r in rc), "every sample matches its deposited column")

    print("summary")
    check(len(tsv(res / "summary.tsv")) == 9, "summary.tsv has 9 samples")

    print(f"\n{'FAILED' if fails else 'OK'}: {len(fails)} failed check(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "results"))
