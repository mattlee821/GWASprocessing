"""Shared helpers for the GWASprocessing test suite.

Usage (from the repo root):
    python3 -m unittest discover -s pipeline/tests -v
pandas, numpy, PyYAML, requests and matplotlib must be importable, and bgzip/tabix must be on PATH for the
end-to-end tests. Liftover tests also need bcftools with the +liftover plugin (GWAS_BCFTOOLS, BCFTOOLS_PLUGINS in .env).
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import unittest
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
PIPELINE = ROOT / "pipeline"
BIN = PIPELINE / "bin"
EXAMPLES = PIPELINE / "tests" / "examples"
SRC = ROOT / "src"

if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

# Tests must never reach the network, even if code under test ignores a test's service override:
# every non-local HTTP(S)/FTP request goes to a dead proxy. Subprocesses inherit this environment.
for _var in ("http_proxy", "https_proxy", "ftp_proxy", "all_proxy", "HTTP_PROXY", "HTTPS_PROXY", "FTP_PROXY", "ALL_PROXY"):
    os.environ[_var] = "http://127.0.0.1:9"
os.environ["no_proxy"] = os.environ["NO_PROXY"] = "127.0.0.1,localhost"
os.environ.setdefault("NXF_OFFLINE", "true")


# Example liftover: on chromosomes 1 and 2, GRCh38 = GRCh37 + 100 for positions 1..1,000,000.
CHAIN_OFFSET = 100


def chain_block(src: str, dst: str, src_size: int, dst_size: int, src_start: int, length: int,
                dst_start: int, chain_id: int, dst_strand: str = "+") -> str:
    """One UCSC chain with a single ungapped block (0-based half-open coordinates)."""
    return (
        f"chain 1000 {src} {src_size} + {src_start} {src_start + length} "
        f"{dst} {dst_size} {dst_strand} {dst_start} {dst_start + length} {chain_id}\n{length}\n\n"
    )


def write_chain(path: Path, blocks: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("".join(blocks))
    return path


def write_offset_chains(genome_dir: Path, offset: int = CHAIN_OFFSET, chroms: tuple[str, ...] = ("1", "2")) -> Path:
    """hg19->hg38 adds `offset`, hg38->hg19 subtracts it, on the first 1 Mb of each chromosome."""
    forward = [chain_block(f"chr{c}", f"chr{c}", 300_000_000, 300_000_000, 0, 1_000_000, offset, i + 1)
               for i, c in enumerate(chroms)]
    backward = [chain_block(f"chr{c}", f"chr{c}", 300_000_000, 300_000_000, offset, 1_000_000, 0, i + 1)
                for i, c in enumerate(chroms)]
    write_chain(genome_dir / "hg19ToHg38.over.chain.gz", forward)
    write_chain(genome_dir / "hg38ToHg19.over.chain.gz", backward)
    return genome_dir


def write_fasta(genome_dir: Path, sequences: dict[str, str], name: str = "GRCh38.fa", width: int = 60) -> Path:
    genome_dir.mkdir(parents=True, exist_ok=True)
    path = genome_dir / name
    with open(path, "w", encoding="utf-8") as handle:
        for chrom, seq in sequences.items():
            handle.write(f">{chrom}\n")
            for start in range(0, len(seq), width):
                handle.write(seq[start:start + width] + "\n")
    return path


def write_grch38_references(root: Path, variants: list[tuple], populations: tuple[str, ...] = ("EUR", "ALL"),
                            paired: list[tuple] | None = None, merges: dict[int, int] | None = None) -> Path:
    """GRCh38 reference arrays as the src/reference_data builders write them.

    variants: (chrom, pos38, rs number, REF, ALT, ALT frequency) -> dbsnp/GRCh38 and 1000genomes/GRCh38 for each
    population; paired: (chrom, pos37, pos38, REF, ALT) -> build_detection sample; merges: {retired rs number: current
    rs number} -> dbsnp/GRCh38/merged_{from,into}.u32."""
    import numpy as np
    import gwas_reference as ref

    db = root / "dbsnp" / "GRCh38"
    db.mkdir(parents=True, exist_ok=True)
    rows = sorted(variants, key=lambda v: (ref.CHROM_CODE[str(v[0])], v[1]))
    for chrom in ref.CHROMS:
        sel = [v for v in rows if str(v[0]) == chrom]
        np.array([v[1] for v in sel], dtype="<u4").tofile(db / f"chr{chrom}.pos.u32")
        np.array(ref.allele_hash([v[3] for v in sel], [v[4] for v in sel]) if sel else [], dtype="<u4").tofile(db / f"chr{chrom}.hash.u32")
        np.array([v[2] for v in sel], dtype="<u4").tofile(db / f"chr{chrom}.rs.u32")
    by_rs = sorted({(v[2], ref.CHROM_CODE[str(v[0])], v[1]) for v in rows if v[2] > 0})
    np.array([r[0] for r in by_rs], dtype="<u4").tofile(db / "rs.u32")
    np.array([r[1] for r in by_rs], dtype="u1").tofile(db / "rs_chrom.u8")
    np.array([r[2] for r in by_rs], dtype="<u4").tofile(db / "rs_pos.u32")
    (db / "meta.json").write_text(json.dumps({"dbsnp_build": "test", "version": 1}), encoding="utf-8")
    if merges:
        retired = sorted(merges)
        np.array(retired, dtype="<u4").tofile(db / "merged_from.u32")
        np.array([merges[r] for r in retired], dtype="<u4").tofile(db / "merged_into.u32")
    kg = root / "1000genomes" / "GRCh38"
    kg.mkdir(parents=True, exist_ok=True)
    snvs = [v for v in rows if len(v[3]) == 1 and len(v[4]) == 1 and 0 < v[5] < 1]
    for pop in populations:
        np.array([ref.CHROM_CODE[str(v[0])] * ref.KEY_SCALE + v[1] for v in snvs], dtype="<i8").tofile(kg / f"{pop}.key.i64")
        np.array([ref.BASES[v[3]] for v in snvs], dtype="u1").tofile(kg / f"{pop}.ref.u1")
        np.array([ref.BASES[v[4]] for v in snvs], dtype="u1").tofile(kg / f"{pop}.alt.u1")
        np.array([v[5] for v in snvs], dtype="<f4").tofile(kg / f"{pop}.af.f4")
    bd = root / "build_detection"  # always written (startup requires it), empty when no pairs are given
    bd.mkdir(parents=True, exist_ok=True)
    (bd / "meta.json").write_text(json.dumps({"snps": len(paired or []), "version": 2}), encoding="utf-8")
    paired = paired or []
    for name, col in (("grch37", 1), ("grch38", 2)):
        order = sorted(paired, key=lambda v: ref.CHROM_CODE[str(v[0])] * ref.KEY_SCALE + v[col])
        np.array([ref.CHROM_CODE[str(v[0])] * ref.KEY_SCALE + v[col] for v in order], dtype="<i8").tofile(bd / f"{name}.key.i64")
        np.array([ref.BASES[v[3]] for v in order], dtype="u1").tofile(bd / f"{name}.ref.u1")
        np.array([ref.BASES[v[4]] for v in order], dtype="u1").tofile(bd / f"{name}.alt.u1")
    return root


def grch38_position(chrom: str, pos19: int, offset: int = CHAIN_OFFSET) -> int:
    """Where write_offset_chains puts a GRCh37 position in GRCh38 (+offset on the first 1 Mb of chromosomes 1 and 2)."""
    return pos19 + offset if str(chrom) in ("1", "2") and 0 < pos19 <= 1_000_000 - offset else pos19


def _chrom(value: Any) -> str:
    text = re.sub(r"^chr", "", str(value), flags=re.I).upper()
    return {"23": "X", "M": "MT"}.get(text, text)


def write_genomes(genome_dir: Path, sites: list[tuple], seed: int = 821, offset: int = CHAIN_OFFSET) -> None:
    """GRCh37 and GRCh38 FASTAs consistent with write_offset_chains, for BCFtools/liftover tests: chromosomes 1 and 2,
    GRCh38 = 100 random bases + GRCh37, with each site's REF base placed at its position in both builds.
    sites: (chrom, pos19, REF)."""
    import numpy as np

    rng = np.random.default_rng(seed)
    length = max([5000] + [int(s[1]) + 1000 for s in sites])
    seq37 = {c: list("".join(rng.choice(list("ACGT"), size=length))) for c in ("1", "2")}
    for chrom, pos19, ref in sites:
        if _chrom(chrom) in seq37 and len(str(ref)) == 1:
            seq37[_chrom(chrom)][int(pos19) - 1] = str(ref).upper()
    pad = "".join(rng.choice(list("ACGT"), size=offset))
    write_fasta(genome_dir, {c: "".join(v) for c, v in seq37.items()}, name="GRCh37.fa")
    write_fasta(genome_dir, {c: pad + "".join(v) for c, v in seq37.items()}, name="GRCh38.fa")


def make_reference(root: Path, lookup_rows: list[tuple] | None = None, population: str = "ALL",
                   offset_chains: bool = True, fasta: dict[str, str] | None = None, genomes: bool = False,
                   genome_sites: list[tuple] | None = None) -> Path:
    """A GRCh38 reference root: dbSNP index, 1000 Genomes frequencies (population, EUR and ALL) and the paired
    build-detection sample, from lookup_rows (CHR, SNP, POS19, A1=ALT, A2=REF, AF1) placed in GRCh38 by the offset
    chains. fasta: GRCh38 sequences to write; genomes: GRCh37 and GRCh38 FASTAs consistent with the chains for
    liftover (REF bases at the lookup sites and at genome_sites, given as (chrom, pos19, REF))."""
    root.mkdir(parents=True, exist_ok=True)
    variants, paired = [], []
    for row in lookup_rows or []:
        chrom, snp, pos19, alt, ref, af = _chrom(row[0]), str(row[1]), int(row[2]), str(row[3]).upper(), str(row[4]).upper(), float(row[5])
        pos38 = grch38_position(chrom, pos19)
        rs = re.fullmatch(r"rs(\d+)", snp, flags=re.I)
        variants.append((chrom, pos38, int(rs.group(1)) if rs else 0, ref, alt, af))
        if len(ref) == 1 and len(alt) == 1:
            paired.append((chrom, pos19, pos38, ref, alt))
    write_grch38_references(root, variants, populations=tuple(dict.fromkeys((population, "EUR", "ALL"))), paired=paired)
    if offset_chains:
        write_offset_chains(root / "genome")
    if genomes:
        write_genomes(root / "genome", [(v[0], int(r[2]), v[3]) for v, r in zip(variants, lookup_rows or [])] + list(genome_sites or []))
    if fasta is not None:
        write_fasta(root / "genome", fasta)
    return root


# The example reference used by the end-to-end tests (formerly a GRCh37 lookup in examples/reference).
EXAMPLE_LOOKUP = [("1", "rs1", 1000, "A", "G", 0.11), ("1", "rs2", 2000, "C", "T", 0.21), ("1", "rs3", 3000, "A", "C", 0.33),
                  ("1", "rs4", 4000, "A", "G", 0.44), ("1", "rs5", 5000, "A", "G", 0.42), ("1", "rs6", 6000, "A", "G", 0.48)]


def bcftools_liftover_available() -> bool:
    tool = os.environ.get("GWAS_BCFTOOLS") or shutil.which("bcftools")
    if not tool:
        return False
    listed = subprocess.run([tool, "plugin", "-l"], capture_output=True, text=True).stdout
    return "liftover" in listed.split()


_EXAMPLE_REFERENCES: dict[bool, Path] = {}


def generic_grch38_table(path: Path) -> Path:
    """examples/tables/generic.tsv with its positions given in GRCh38 (POS38 column), so no liftover is needed."""
    return write_table(path, ["CHR", "POS38", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele", "BETA", "SE", "P", "N"],
                       [[1, 1100, "rs1", "A", "G", 0.10, 0.20, 0.05, 0.001, 1000],
                        [1, 2100, "rs2", "C", "T", 0.20, -0.10, 0.04, 0.5, 1000]])


needs_liftover = unittest.skipUnless(bcftools_liftover_available(), "bcftools with the +liftover plugin is required "
                                     "(GWAS_BCFTOOLS, BCFTOOLS_PLUGINS in .env)")


def example_reference(genomes: bool = False) -> Path:
    """The example reference (EXAMPLE_LOOKUP), built once per test process; genomes=True adds GRCh37 and GRCh38 FASTAs
    consistent with the offset chains, for inputs that need liftover."""
    if genomes not in _EXAMPLE_REFERENCES:
        import tempfile

        root = Path(tempfile.mkdtemp(prefix="gwas_example_reference_"))
        _EXAMPLE_REFERENCES[genomes] = make_reference(root, lookup_rows=EXAMPLE_LOOKUP, genomes=genomes)
    return _EXAMPLE_REFERENCES[genomes]


def base_record(raw: Path, output_prefix: Path, **overrides: Any) -> dict[str, Any]:
    record = {
        "GWAS_location": str(raw), "phenotype": "BMI", "ancestry": "ALL", "sex": "combined",
        "author": "Test", "year": "2026", "PMID": "1", "study_id": "test_2026_1", "output_id": "test",
        "row_hash": "abc", "source_type": "local", "raw_path": str(raw), "output_prefix": str(output_prefix),
        "standardise_dir": str(output_prefix.parent), "metadata": {}, "input_build": "", "population": "",
        "delimiter": "", "format": "", "study_yaml": "",
    }
    record.update(overrides)
    return record


def run_script(script: str, *args: Any, env: dict[str, str] | None = None, check: bool = True,
               cwd: Path | None = None) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(BIN / script), *map(str, args)]
    proc = subprocess.run(cmd, cwd=cwd or ROOT, text=True, capture_output=True, env=env)
    if check and proc.returncode != 0:
        raise AssertionError(f"{script} exited {proc.returncode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc


@dataclass
class StandardiseResult:
    returncode: int
    stderr: str
    record: dict[str, Any] = field(default_factory=dict)
    qc: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    text: str = ""

    @property
    def header(self) -> list[str]:
        return [line for line in self.text.splitlines() if line.startswith("#")]

    @property
    def rows(self) -> list[dict[str, Any]]:
        return parse_vcf_rows(self.text)

    def by_id(self) -> dict[str, dict[str, Any]]:
        return {row["ID"]: row for row in self.rows}

    def by_snpid(self) -> dict[str, dict[str, Any]]:
        """Rows by INFO/SNPID, the study's own identifier (ID is set only from dbSNP)."""
        return {row["INFO"].get("SNPID", ""): row for row in self.rows}

    @property
    def sample(self) -> dict[str, str]:
        """The ##SAMPLE header line as a dict (quoted values unquoted)."""
        return parse_sample_line(self.header)


def parse_sample_line(header: list[str]) -> dict[str, str]:
    import re as _re
    line = next((h for h in header if h.startswith("##SAMPLE=<")), "")
    pairs = _re.findall(r'([A-Za-z]+)=("(?:[^"\\]|\\.)*"|[^,>]*)', line[len("##SAMPLE=<"):])
    out = {}
    for key, value in pairs:
        if value.startswith('"'):
            value = _re.sub(r'\\(.)', r'\1', value[1:-1])
        out[key] = value
    return out


def parse_vcf_rows(text: str) -> list[dict[str, Any]]:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        info = {}
        for item in parts[7].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info[key] = value
        fmt = dict(zip(parts[8].split(":"), parts[9].split(":"))) if len(parts) > 9 else {}
        rows.append({"CHROM": parts[0], "POS": int(parts[1]), "ID": parts[2], "REF": parts[3], "ALT": parts[4],
                     "INFO": info, "FORMAT": fmt, "line": line})
    return rows


def read_vcf_text(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return handle.read()


def run_standardise(tmp: Path, raw: Path, reference_root: Path, *,
                    config_root: Path | None = None, study_yaml: str = "", force: bool = False,
                    env: dict[str, str] | None = None, name: str = "case", **record_overrides: Any) -> StandardiseResult:
    case = tmp / name
    case.mkdir(parents=True, exist_ok=True)
    output_prefix = case / "standardise" / "standardised"
    stage = case / "stage.json"
    stage.write_text(json.dumps(base_record(raw, output_prefix, **record_overrides)), encoding="utf-8")
    out = case / "record.json"
    args = ["--stage-record", stage, "--config-root", config_root or (case / "no_manifests"),
            "--reference-root", reference_root,
            "--force", "true" if force else "false", "--output", out]
    if study_yaml:
        args += ["--study-yaml", study_yaml]
    proc = run_script("standardise_gwas.py", *args, env=env, check=False)
    result = StandardiseResult(returncode=proc.returncode, stderr=proc.stderr)
    if proc.returncode == 0:
        result.record = json.loads(out.read_text(encoding="utf-8"))
        result.qc = json.loads(Path(result.record["qc_json"]).read_text(encoding="utf-8"))
        result.metadata = json.loads(Path(result.record["metadata_json"]).read_text(encoding="utf-8"))
        result.text = read_vcf_text(Path(result.record["standard_vcf"]))
    return result


def write_table(path: Path, header: list[str], rows: list[list[Any]], sep: str = "\t") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [sep.join(header)] + [sep.join("" if v is None else str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class LocalHttpServer:
    """Serve a directory over HTTP on 127.0.0.1 for download tests (context manager)."""

    def __init__(self, directory: Path):
        import functools
        import http.server
        import threading

        class QuietHandler(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args, **kwargs):
                pass

        handler = functools.partial(QuietHandler, directory=str(directory))
        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.server.server_address[1]}"

    def __enter__(self) -> "LocalHttpServer":
        self.thread.start()
        return self

    def __exit__(self, *exc: Any) -> None:
        self.server.shutdown()
        self.server.server_close()


def write_stub(bin_dir: Path, name: str, body: str) -> Path:
    """Create an executable shell stub (e.g. a fake `aws` or `nextflow`) in bin_dir."""
    bin_dir.mkdir(parents=True, exist_ok=True)
    path = bin_dir / name
    path.write_text("#!/usr/bin/env bash\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path
