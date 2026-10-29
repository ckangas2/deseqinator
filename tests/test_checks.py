"""
Unit tests: one planted failure per check, plus the passing case.

    bash tests/run_tests.sh        (or: python -m pytest tests -q)

No STAR, R or network. Each script is run as a subprocess, exactly as the
pipeline runs it, and judged on exit code and output files.
"""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
DATA = ROOT / "tests" / "data"


def run(script, *args, cwd=None):
    p = subprocess.run([sys.executable, "-I", str(BIN / script), *map(str, args)],
                       capture_output=True, text=True, cwd=cwd)
    return p.returncode, p.stdout + p.stderr


def write(path, text):
    path.write_text(text)
    return path


# ------------------------------------------------------------------ design

SHEET_HDR = "sample,fastq_1,fastq_2,cell_line,condition,replicate,batch\n"


def sheet(tmp_path, rows, name="ss.csv"):
    return write(tmp_path / name, SHEET_HDR + "".join(",".join(r) + "\n" for r in rows))


def base_rows(batch=lambda cond: "b1"):
    rows = []
    for cond in ("Mock", "Virus"):
        for rep in "123":
            rows.append([f"A_{cond}_{rep}", f"A_{cond}_{rep}.fq.gz", "", "A", cond, rep, batch(cond)])
    return rows


def contrasts(tmp_path, body):
    return write(tmp_path / "ct.csv", "id,cell_line,reference,target\n" + body)


def test_design_passes_and_scopes_to_one_line(tmp_path):
    rows = base_rows() + [["B_Mock_1", "b.fq.gz", "", "B", "Mock", "1", "b9"],
                          ["A_Rux_1", "r.fq.gz", "", "A", "Rux", "1", "b1"]]
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "c1,A,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 0, out
    v = list(csv.DictReader(open(tmp_path / "d/contrasts.validated.tsv"), delimiter="\t"))
    assert v[0]["design"] == "~ condition"          # batch constant within contrast: dropped
    assert set(v[0]["samples"].split(",")) == {r[0] for r in base_rows()}
    assert "Rux" in out and "B_Mock_1" in out       # exclusions are reported


def test_design_rejects_underpowered(tmp_path):
    rows = base_rows()[:2] + base_rows()[3:]
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "c1,A,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 1 and "need >= 3 per arm" in out


def test_design_rejects_confounded_batch(tmp_path):
    rows = base_rows(batch=lambda c: f"run_{c}")
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "c1,A,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 1 and "confounded" in out


def test_design_keeps_partially_crossed_batch(tmp_path):
    rows = base_rows(batch=lambda c: "b1")
    rows[0][6] = rows[3][6] = "b2"          # one mock and one infected in a second batch
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "c1,A,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 0, out
    v = next(csv.DictReader(open(tmp_path / "d/contrasts.validated.tsv"), delimiter="\t"))
    assert v["design"] == "~ batch + condition"


def test_design_star_contrast_adds_cell_line(tmp_path):
    rows = base_rows() + [[f"B_{c}_{r}", "x.fq.gz", "", "B", c, r, "b1"]
                          for c in ("Mock", "Virus") for r in "123"]
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "all,*,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 0, out
    v = next(csv.DictReader(open(tmp_path / "d/contrasts.validated.tsv"), delimiter="\t"))
    assert v["design"] == "~ cell_line + condition"


def test_design_rejects_line_confounded_with_condition(tmp_path):
    rows = [r for r in base_rows() if r[4] == "Mock"] + \
           [[f"B_Virus_{r}", "x.fq.gz", "", "B", "Virus", r, "b1"] for r in "123"]
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "all,*,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 1 and "confounded" in out


@pytest.mark.parametrize("mutate,needle", [
    (lambda r: r[1].__setitem__(0, r[0][0]), "duplicate sample name"),
    (lambda r: r[1].__setitem__(5, "1"), "share cell_line/condition/replicate"),
    (lambda r: r[0].__setitem__(2, "x_R2.fq.gz"), "mixed single-end"),
])
def test_design_rejects_sheet_errors(tmp_path, mutate, needle):
    rows = base_rows()
    mutate(rows)
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, rows),
                  "--contrasts", contrasts(tmp_path, "c1,A,Mock,Virus\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 1 and needle in out


def test_design_rejects_unknown_condition(tmp_path):
    rc, out = run("check_design.py", "--samplesheet", sheet(tmp_path, base_rows()),
                  "--contrasts", contrasts(tmp_path, "c1,A,Mock,SARS-CoV-2\n"),
                  "--out-dir", tmp_path / "d")
    assert rc == 1 and "matched exactly" in out


# ------------------------------------------------------------------ strandedness

def rpg(tmp_path, name, col3, col4, genes=50):
    lines = ["N_unmapped\t10\t10\t10", "N_multimapping\t5\t5\t5",
             f"N_noFeature\t1\t{col4 * genes}\t{col3 * genes}", "N_ambiguous\t0\t0\t0"]
    lines += [f"G{i}\t{col3 + col4}\t{col3}\t{col4}" for i in range(genes)]
    return write(tmp_path / f"{name}.ReadsPerGene.out.tab", "\n".join(lines) + "\n")


@pytest.mark.parametrize("c3,c4,want,code", [
    (95, 5, "forward", 1), (4, 96, "reverse", 2), (50, 52, "unstranded", 0)])
def test_strand_calls(tmp_path, c3, c4, want, code):
    files = [rpg(tmp_path, f"s{i}", c3, c4) for i in range(3)]
    rc, out = run("infer_strandedness.py", *files, "--out-tsv", tmp_path / "s.tsv",
                  "--out-call", tmp_path / "s.txt")
    assert rc == 0, out
    assert (tmp_path / "s.txt").read_text().strip() == want
    assert f"-s {code}" in out


def test_strand_ambiguous_fails(tmp_path):
    files = [rpg(tmp_path, "s1", 70, 30)]
    rc, out = run("infer_strandedness.py", *files, "--out-tsv", tmp_path / "s.tsv",
                  "--out-call", tmp_path / "s.txt")
    assert rc == 1 and "no clear strandedness" in out
    assert not (tmp_path / "s.txt").exists()


def test_strand_disagreement_fails(tmp_path):
    files = [rpg(tmp_path, "s1", 95, 5), rpg(tmp_path, "s2", 5, 95)]
    rc, out = run("infer_strandedness.py", *files, "--out-tsv", tmp_path / "s.tsv",
                  "--out-call", tmp_path / "s.txt")
    assert rc == 1 and "disagree" in out


def test_strand_declared_mismatch_fails(tmp_path):
    files = [rpg(tmp_path, "s1", 5, 95)]
    rc, out = run("infer_strandedness.py", *files, "--declared", "forward",
                  "--out-tsv", tmp_path / "s.tsv", "--out-call", tmp_path / "s.txt")
    assert rc == 1 and "--strandedness forward was given" in out


def test_strand_too_few_reads_fails(tmp_path):
    files = [rpg(tmp_path, "s1", 1, 9, genes=5)]
    rc, out = run("infer_strandedness.py", *files, "--out-tsv", tmp_path / "s.tsv",
                  "--out-call", tmp_path / "s.txt")
    assert rc == 1


# ------------------------------------------------------------------ featureCounts

def fc(tmp_path, assigned, nofeat):
    hdr = "Geneid\tChr\tStart\tEnd\tStrand\tLength\tbam/a.bam\tbam/b.bam\n"
    write(tmp_path / "fc.txt", "# Program:featureCounts\n" + hdr +
          "G1\tchr1\t1\t10\t+\t10\t5\t6\nG2\tchr1\t20\t30\t+\t11\t7\t8\n")
    write(tmp_path / "fc.txt.summary",
          "Status\tbam/a.bam\tbam/b.bam\n"
          f"Assigned\t{assigned[0]}\t{assigned[1]}\n"
          f"Unassigned_NoFeatures\t{nofeat[0]}\t{nofeat[1]}\n"
          "Unassigned_Ambiguity\t1\t1\n")
    return tmp_path / "fc.txt"


def test_counts_tidy_and_pass(tmp_path):
    f = fc(tmp_path, (80, 75), (20, 25))
    rc, out = run("tidy_counts.py", "--featurecounts", f, "--out-matrix", tmp_path / "m.tsv",
                  "--out-qc", tmp_path / "q.tsv")
    assert rc == 0, out
    assert (tmp_path / "m.tsv").read_text().splitlines()[0] == "gene_id\ta\tb"


def test_counts_low_assignment_fails(tmp_path):
    f = fc(tmp_path, (80, 8), (20, 92))
    rc, out = run("tidy_counts.py", "--featurecounts", f, "--out-matrix", tmp_path / "m.tsv",
                  "--out-qc", tmp_path / "q.tsv")
    assert rc == 1 and "b:" in out and "Unassigned_NoFeatures" in out


# ------------------------------------------------------------------ viral gate

def loads(tmp_path, vals):
    for s, (t, v) in vals.items():
        write(tmp_path / f"{s}.viral.tsv", f"{s}\t{t}\t{v}\n")
    return [tmp_path / f"{s}.viral.tsv" for s in vals]


GATE_SHEET = "sample,condition\nm1,Mock\nm2,Mock\ni1,Virus\ni2,Virus\n"


def test_gate_passes(tmp_path):
    f = loads(tmp_path, {"m1": (1000, 0), "m2": (1000, 0), "i1": (1000, 200), "i2": (1000, 150)})
    rc, out = run("viral_gate.py", *f, "--samplesheet", write(tmp_path / "s.csv", GATE_SHEET),
                  "--out", tmp_path / "v.tsv")
    assert rc == 0, out


def test_gate_viral_mock_fails(tmp_path):
    f = loads(tmp_path, {"m1": (1000, 200), "m2": (1000, 0), "i1": (1000, 0), "i2": (1000, 150)})
    rc, out = run("viral_gate.py", *f, "--samplesheet", write(tmp_path / "s.csv", GATE_SHEET),
                  "--out", tmp_path / "v.tsv")
    assert rc == 1 and "m1 is labelled Mock" in out and "i1 is labelled Virus" in out


def test_gate_missing_sample_fails(tmp_path):
    f = loads(tmp_path, {"m1": (1000, 0), "m2": (1000, 0), "i1": (1000, 200)})
    rc, out = run("viral_gate.py", *f, "--samplesheet", write(tmp_path / "s.csv", GATE_SHEET),
                  "--out", tmp_path / "v.tsv")
    assert rc == 1 and "no viral load reported" in out


# ------------------------------------------------------------------ reference

def test_reference_passes_on_test_data(tmp_path):
    contigs = write(tmp_path / "c.txt", ">chr1 x\n>chr2 x\n")
    rc, out = run("validate_reference.py", "--host-contigs", contigs,
                  "--host-gtf", DATA / "host.gtf.gz", "--viral-fasta", DATA / "virus.fa.gz")
    assert rc == 0, out


def test_reference_contig_style_mismatch_fails(tmp_path):
    contigs = write(tmp_path / "c.txt", ">chr1\n>chr2\n")
    rc, out = run("validate_reference.py", "--host-contigs", contigs,
                  "--host-gtf", DATA / "host_nochr.gtf.gz", "--viral-fasta", DATA / "virus.fa.gz")
    assert rc == 1 and "Ensembl annotation on a UCSC/GENCODE genome" in out


def test_reference_viral_name_collision_fails(tmp_path):
    contigs = write(tmp_path / "c.txt", ">chr1\n>chr2\n>MN985325.1\n")
    rc, out = run("validate_reference.py", "--host-contigs", contigs,
                  "--host-gtf", DATA / "host.gtf.gz", "--viral-fasta", DATA / "virus.fa.gz")
    assert rc == 1 and "also used by the host" in out


# ------------------------------------------------------------------ reconcile

def our_counts(tmp_path):
    """Our matrix = the deposited truth, keyed by gene_id, named by sample."""
    rows = list(csv.reader(open(DATA / "deposited_counts.tsv"), delimiter="\t"))
    ss = list(csv.DictReader(open(DATA / "samplesheet.csv")))
    by_dep = {r["deposited_id"]: r["sample"] for r in ss}
    ann = {}
    with gzip.open(DATA / "host.gtf.gz", "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[2] == "gene":
                gid = f[8].split('"')[1]; name = f[8].split('"')[3]
                ann[name] = gid
    write(tmp_path / "ann.tsv", "gene_id\tgene_id_unversioned\tgene_name\tgene_type\tchrom\n" +
          "".join(f"{g}\t{g.split('.')[0]}\t{n}\tprotein_coding\tchr1\n" for n, g in ann.items()))
    hdr = ["gene_id"] + [by_dep[c] for c in rows[0][1:]]
    write(tmp_path / "ours.tsv", "\t".join(hdr) + "\n" +
          "".join("\t".join([ann[r[0]]] + r[1:]) + "\n" for r in rows[1:]))
    return tmp_path / "ours.tsv", tmp_path / "ann.tsv"


def test_reconcile_matches(tmp_path):
    ours, ann = our_counts(tmp_path)
    rc, out = run("reconcile_counts.py", "--ours", ours, "--annotation", ann,
                  "--deposited", DATA / "deposited_counts.tsv",
                  "--samplesheet", DATA / "samplesheet.csv", "--min-rho", "0.8",
                  "--out", tmp_path / "r.tsv", "--out-matrix", tmp_path / "rm.tsv")
    assert rc == 0, out


def test_reconcile_detects_swap(tmp_path):
    ours, ann = our_counts(tmp_path)
    rc, out = run("reconcile_counts.py", "--ours", ours, "--annotation", ann,
                  "--deposited", DATA / "deposited_counts.tsv",
                  "--samplesheet", DATA / "samplesheet_depswap.csv",
                  "--out", tmp_path / "r.tsv", "--out-matrix", tmp_path / "rm.tsv")
    assert rc == 1 and "Sample swap?" in out


# ------------------------------------------------------------------ gene table

def test_gene_table(tmp_path):
    rc, out = run("gtf_gene_table.py", "--gtf", DATA / "host.gtf.gz", "--out", tmp_path / "g.tsv")
    assert rc == 0, out
    rows = list(csv.DictReader(open(tmp_path / "g.tsv"), delimiter="\t"))
    assert len(rows) == 200 and rows[0]["gene_name"] == "GENE001"
    assert rows[0]["gene_id"] != rows[0]["gene_id_unversioned"]


# ------------------------------------------------------------------ test data

def test_fastq_records_are_well_formed():
    """Every record: @header, sequence, '+', quality of the same length.

    fastp rejects a mismatched record outright; a generator bug here once
    produced 74 bp reads with 75 bp quality strings.
    """
    for fq in sorted(DATA.glob("*.fastq.gz")):
        with gzip.open(fq, "rt") as fh:
            lines = fh.read().splitlines()
        assert len(lines) % 4 == 0, fq.name
        for i in range(0, len(lines), 4):
            h, s, plus, q = lines[i:i + 4]
            assert h.startswith("@") and plus == "+", f"{fq.name} record {i // 4}"
            assert len(s) == len(q) and s, f"{fq.name} {h}: seq {len(s)} vs qual {len(q)}"
