#!/usr/bin/env python3
"""
Deterministically regenerate tests/data/: a toy host + virus reference,
nine single-end reverse-stranded libraries with planted differential
expression, gene sets built on the planted genes, a "deposited" count
matrix, and the broken variants CI uses for negative cases.

    python3 tests/make_test_data.py tests/data

Truth is written to truth.tsv; tests/check_test_results.py reads it.
"""
import gzip
import random
import sys
from pathlib import Path

import numpy as np

SEED = 20201
READ_LEN = 75
N_GENES = 200
EXON, INTRON, SPACER = 300, 400, 500
VIRAL_LEN = 5000
VIRAL_CONTIG = "MN985325.1"
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def gz_write(path, text):
    # mtime=0 so regenerating gives byte-identical files
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write(text.encode())


def wrap(seq, n=60):
    return "\n".join(seq[i:i + n] for i in range(0, len(seq), n))


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rnd = random.Random(SEED)

    # ------------------------------------------------------------ genome
    genes = []          # (gene_id, name, chrom, strand, exon1, exon2)
    chroms = {}
    per_chrom = N_GENES // 2
    gi = 0
    for chrom in ("chr1", "chr2"):
        pos = SPACER
        for _ in range(per_chrom):
            gi += 1
            strand = "+" if gi % 2 else "-"
            e1 = (pos + 1, pos + EXON)
            e2 = (pos + EXON + INTRON + 1, pos + 2 * EXON + INTRON)
            genes.append((f"ENSGT{gi:011d}.{1 + gi % 3}", f"GENE{gi:03d}", chrom, strand, e1, e2))
            pos += 2 * EXON + INTRON + SPACER
        seq = "".join(rnd.choice("ACGT") for _ in range(pos + SPACER))
        # canonical splice sites so STAR's annotated junctions look real
        s = list(seq)
        for g in genes:
            if g[2] != chrom:
                continue
            e1, e2 = g[4], g[5]
            donor, acceptor = e1[1], e2[0] - 1     # 1-based last exon base / last intron base
            if g[3] == "+":
                s[donor:donor + 2] = "GT"; s[acceptor - 2:acceptor] = "AG"
            else:
                s[donor:donor + 2] = "CT"; s[acceptor - 2:acceptor] = "AC"
        chroms[chrom] = "".join(s)
    virus = "".join(rnd.choice("ACGT") for _ in range(VIRAL_LEN))

    fa = "".join(f">{c} synthetic host\n{wrap(s)}\n" for c, s in chroms.items())
    gz_write(out / "host.fa.gz", fa)
    gz_write(out / "virus.fa.gz", f">{VIRAL_CONTIG} synthetic virus\n{wrap(virus)}\n")

    def gtf_text(rename=None):
        rows = []
        for gid, name, chrom, strand, e1, e2 in genes:
            c = rename(chrom) if rename else chrom
            attrs = f'gene_id "{gid}"; gene_name "{name}"; gene_type "protein_coding";'
            tid = gid.replace("ENSGT", "ENSTT")
            rows.append(f"{c}\tTEST\tgene\t{e1[0]}\t{e2[1]}\t.\t{strand}\t.\t{attrs}")
            rows.append(f'{c}\tTEST\ttranscript\t{e1[0]}\t{e2[1]}\t.\t{strand}\t.\t{attrs} transcript_id "{tid}";')
            for n, (a, b) in enumerate((e1, e2), 1):
                rows.append(f'{c}\tTEST\texon\t{a}\t{b}\t.\t{strand}\t.\t{attrs} transcript_id "{tid}"; exon_number {n};')
        return "\n".join(rows) + "\n"

    gz_write(out / "host.gtf.gz", gtf_text())
    # Ensembl-style contig names against a chr-style FASTA: must fail validation
    gz_write(out / "host_nochr.gtf.gz", gtf_text(lambda c: c.replace("chr", "")))

    def transcript(g):
        _, _, chrom, strand, e1, e2 = g
        s = chroms[chrom]
        t = s[e1[0] - 1:e1[1]] + s[e2[0] - 1:e2[1]]
        return t if strand == "+" else revcomp(t)

    tx = [transcript(g) for g in genes]

    # ------------------------------------------------------------ design
    samples = []
    for rep in (1, 2, 3):
        samples.append(dict(sample=f"lineA_mock_{rep}", cell_line="lineA", condition="Mock",
                            replicate=rep, batch="b1", deposited_id=f"Series1_lineA_Mock_{rep}"))
    for rep in (1, 2, 3):
        samples.append(dict(sample=f"lineA_inf_{rep}", cell_line="lineA", condition="SARS-CoV-2",
                            replicate=rep, batch="b1", deposited_id=f"Series1_lineA_SARS-CoV-2_{rep}"))
    samples.append(dict(sample="lineA_rux_1", cell_line="lineA", condition="SARS-CoV-2_Rux",
                        replicate=1, batch="b2", deposited_id="Series2_lineA_SARS-CoV-2_Rux_1"))
    samples.append(dict(sample="lineB_mock_1", cell_line="lineB", condition="Mock",
                        replicate=1, batch="b3", deposited_id="Series3_lineB_Mock_1"))
    samples.append(dict(sample="lineB_inf_1", cell_line="lineB", condition="SARS-CoV-2",
                        replicate=1, batch="b3", deposited_id="Series3_lineB_SARS-CoV-2_1"))

    # ------------------------------------------------------------ expression
    base = rng.lognormal(mean=5.0, sigma=0.9, size=N_GENES)
    base = base / base.sum() * 24000
    eligible = [i for i in range(N_GENES) if base[i] >= 60]
    picks = rng.choice(eligible, size=40, replace=False)
    up, down = sorted(picks[:20].tolist()), sorted(picks[20:].tolist())
    null_pool = [i for i in range(N_GENES) if i not in set(up) | set(down) and base[i] >= 40]
    null_pick = rng.choice(null_pool, size=30, replace=False).tolist()
    # lineB differs strongly at baseline: pooling it into lineA's contrast would show
    lineB_shift = rng.choice(N_GENES, size=50, replace=False)

    dispersion = 0.03
    truth_counts = {}
    for smp in samples:
        mu = base.copy()
        if smp["cell_line"] == "lineB":
            mu[lineB_shift] *= 6
        if smp["condition"] != "Mock":
            mu[up] *= 4
            mu[down] /= 4
        mu *= rng.uniform(0.8, 1.25)              # library size
        r = 1 / dispersion
        counts = rng.negative_binomial(r, r / (r + mu))
        truth_counts[smp["sample"]] = counts
        viral_n = int(round(counts.sum() * 0.25)) if smp["condition"] != "Mock" else 0

        lines = []
        k = 0
        def emit(seq):
            nonlocal k
            k += 1
            if rnd.random() < 0.05:                # short insert, adapter read-through
                cut = rnd.randint(40, READ_LEN - 15)
                seq = (seq[:cut] + ADAPTER)[:READ_LEN]
            lines.append(f"@{smp['sample']}.{k}\n{seq}\n+\n{'I' * len(seq)}\n")
        for g, n in enumerate(counts):
            t = tx[g]
            for _ in range(n):
                p = rnd.randint(0, len(t) - READ_LEN)
                emit(revcomp(t[p:p + READ_LEN]))        # reverse-stranded (dUTP)
        for _ in range(viral_n):
            p = rnd.randint(0, VIRAL_LEN - READ_LEN)
            emit(revcomp(virus[p:p + READ_LEN]))
        rnd.shuffle(lines)
        gz_write(out / f"{smp['sample']}.fastq.gz", "".join(lines))
        smp["fastq_1"] = f"{smp['sample']}.fastq.gz"
        smp["fastq_2"] = ""

    # ------------------------------------------------------------ sheets
    cols = ["sample", "fastq_1", "fastq_2", "cell_line", "condition", "replicate", "batch",
            "deposited_id"]
    def sheet(rows):
        return ",".join(cols) + "\n" + "".join(
            ",".join(str(r[c]) for c in cols) + "\n" for r in rows)

    (out / "samplesheet.csv").write_text(sheet(samples))
    (out / "contrasts.csv").write_text(
        "id,cell_line,reference,target\nlineA_inf_vs_mock,lineA,Mock,SARS-CoV-2\n")

    # negative cases for CI
    (out / "contrasts_underpowered.csv").write_text(
        "id,cell_line,reference,target\nlineB_inf_vs_mock,lineB,Mock,SARS-CoV-2\n")
    conf = [dict(s) for s in samples]
    for s in conf:
        if s["cell_line"] == "lineA" and s["condition"] in ("Mock", "SARS-CoV-2"):
            s["batch"] = "batch_" + s["condition"]
    (out / "samplesheet_confounded.csv").write_text(sheet(conf))
    mis = [dict(s) for s in samples]
    a = next(s for s in mis if s["sample"] == "lineA_mock_1")
    b = next(s for s in mis if s["sample"] == "lineA_inf_1")
    a["fastq_1"], b["fastq_1"] = b["fastq_1"], a["fastq_1"]   # files swapped, labels kept
    (out / "samplesheet_mislabelled.csv").write_text(sheet(mis))
    swp = [dict(s) for s in samples]
    a = next(s for s in swp if s["sample"] == "lineA_mock_2")
    b = next(s for s in swp if s["sample"] == "lineA_inf_2")
    a["deposited_id"], b["deposited_id"] = b["deposited_id"], a["deposited_id"]
    (out / "samplesheet_depswap.csv").write_text(sheet(swp))

    # ------------------------------------------------------------ gene sets
    name = lambda i: genes[i][1]
    gid = lambda i: genes[i][0].split(".")[0]
    up_set = [name(i) for i in up[:15]] + [name(i) for i in null_pick[:5]]
    down_set = [name(i) for i in down[:15]] + [name(i) for i in null_pick[5:10]]
    null_set = [name(i) for i in null_pick[10:30]]
    sets = {"HALLMARK_TEST_UP": up_set, "HALLMARK_TEST_DOWN": down_set,
            "HALLMARK_TEST_NULL": null_set}
    (out / "test_hallmarks.gmt").write_text("".join(
        f"{k}\tsynthetic\t" + "\t".join(v) + "\n" for k, v in sets.items()))
    # same sets keyed by Ensembl ID: must fail the coverage check
    by_name = {genes[i][1]: gid(i) for i in range(N_GENES)}
    (out / "test_hallmarks_ensembl.gmt").write_text("".join(
        f"{k}\tsynthetic\t" + "\t".join(by_name[g] for g in v) + "\n" for k, v in sets.items()))

    # ------------------------------------------------------------ deposited + truth
    dep_cols = [s["deposited_id"] for s in samples]
    with open(out / "deposited_counts.tsv", "w") as fh:
        fh.write("gene\t" + "\t".join(dep_cols) + "\n")
        for i in range(N_GENES):
            fh.write(name(i) + "\t" + "\t".join(
                str(truth_counts[s["sample"]][i]) for s in samples) + "\n")
    with open(out / "truth.tsv", "w") as fh:
        fh.write("gene_id\tgene_name\tplanted\tbase_mean\n")
        for i, g in enumerate(genes):
            p = "up" if i in up else "down" if i in down else "none"
            fh.write(f"{g[0]}\t{g[1]}\t{p}\t{base[i]:.1f}\n")

    total = sum(int(c.sum()) for c in truth_counts.values())
    print(f"wrote {out}: {N_GENES} genes, {len(samples)} libraries, {total:,} host reads, "
          f"{len(up)} up / {len(down)} down planted")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
