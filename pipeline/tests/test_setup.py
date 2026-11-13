"""Setup and reference data: src/source.sh (.env, tools, container) and src/reference_data/ (download, builders,
validator), and the reference identity that decides when references have changed.

1 test; each README rule is a subTest. Scripts run inside a temporary copy of the repository with stub executables
(micromamba, apptainer, docker, wget) recording how they were called; network tools and installers are replaced by
stubs that fail. Synthetic data only.
"""

from __future__ import annotations

import bz2
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import (ROOT, SRC, bcftools_liftover_available, write_chain, write_genomes, write_grch38_references,
                     write_offset_chains, write_stub)

import gwas_alleles as alleles
import gwas_pipeline_common as common

BASH = shutil.which("bash")
BLOCKED_TOOLS = ("wget", "curl", "brew", "apt-get", "yum", "sudo", "conda", "mamba", "aws", "gsutil")
POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")
REFDATA = SRC / "reference_data"
PATH_KEYS = ("GWAS_WORK_ROOT", "GWAS_REFERENCE_ROOT", "GWAS_TEMP_ROOT", "GWAS_LOG_ROOT", "GWAS_NEXTFLOW_WORK",
             "GWAS_NEXTFLOW_LOG")
PATH_DEFAULTS = ("GWAS", "references", "temp", "temp/logs", "temp/nextflow_work", "temp/.nextflow.log")


def sandbox_env(tmp: Path, *bin_dirs: Path, **extra: str) -> dict:
    """PATH with the given stub folders first, then stubs that block network and installers, bash, and python3 only
    from the test's interpreter; GWAS_* and related settings from the caller's environment are removed."""
    guard = tmp / "guard_bin"
    for tool in BLOCKED_TOOLS:
        if not (guard / tool).exists():
            write_stub(guard, tool, f'echo "BLOCKED {tool} $*" >&2\nexit 97\n')
    bash_dir = tmp / "bash_bin"
    bash_dir.mkdir(exist_ok=True)
    if not (bash_dir / "bash").exists():
        (bash_dir / "bash").symlink_to(BASH)
    python_dir = tmp / "python_bin"
    python_dir.mkdir(exist_ok=True)
    if not (python_dir / "python3").exists():
        write_stub(python_dir, "python3", f'exec "{sys.executable}" "$@"\n')
    path = [*map(str, bin_dirs), str(guard), str(bash_dir), str(python_dir), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("GWAS_", "SLURM_", "NEXTFLOW", "GWASPROCESSING_", "PLINK2", "BCFTOOLS", "MAMBA",
                                "OPENGWAS_"))}
    env.update({"PATH": os.pathsep.join(path), "PIP_NO_INDEX": "1", "PIP_DISABLE_PIP_VERSION_CHECK": "1"})
    env.update(extra)
    return env


def fake_repo(parent: Path, dotenv: str | None = None) -> Path:
    """A copy of the parts of the repository setup uses."""
    repo = parent / "repo"
    (repo / "src").mkdir(parents=True)
    (repo / "pipeline").mkdir()
    (repo / "pipeline" / "main.nf").write_text("// placeholder\n", encoding="utf-8")
    shutil.copytree(ROOT / "pipeline" / "bin", repo / "pipeline" / "bin", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy(SRC / "source.sh", repo / "src" / "source.sh")
    shutil.copytree(SRC / "reference_data", repo / "src" / "reference_data")
    shutil.copytree(SRC / "lib", repo / "src" / "lib")
    shutil.copy(ROOT / "environment.yml", repo / "environment.yml")
    shutil.copy(ROOT / ".env.example", repo / ".env.example")
    if dotenv is not None:
        (repo / ".env").write_text(dotenv, encoding="utf-8")
    return repo


def settings(repo: Path) -> dict[str, str]:
    lines = (repo / ".env").read_text(encoding="utf-8").splitlines()
    return dict(line.split("=", 1) for line in lines if "=" in line and not line.lstrip().startswith("#"))


def run(script: Path, *args) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True, timeout=300)


class Setup:
    """src/source.sh in a fresh fake repository."""

    def __init__(self, tmp: Path, dotenv: str | None = None, parent: Path | None = None):
        self.tmp = tmp
        tmp.mkdir(parents=True, exist_ok=True)
        self.repo = fake_repo(parent or tmp, dotenv)
        self.mamba_log = tmp / "mamba.log"

    def source(self, *args: str, bin_dirs: tuple = (), **extra: str) -> subprocess.CompletedProcess:
        return subprocess.run([BASH, str(self.repo / "src" / "source.sh"), *args], capture_output=True, text=True,
                              env=sandbox_env(self.tmp, *bin_dirs, **extra), timeout=120)

    def tools(self, *args: str, bin_dirs: tuple = ()) -> subprocess.CompletedProcess:
        return self.source("--tools-only", *args, bin_dirs=bin_dirs, GWAS_DATA_HOME=str(self.tmp / "data"))

    def install_stubs(self) -> Path:
        """A stub micromamba that records its calls and creates an environment of stub tools (its bcftools lists the
        liftover plugin, so nothing is built)."""
        write_stub(self.repo / ".tools" / "bin", "micromamba", f'''echo "micromamba $*" >> "{self.mamba_log}"
prefix=""; while [[ $# -gt 0 ]]; do [[ "$1" == "-p" ]] && prefix="$2"; shift; done
mkdir -p "$prefix/conda-meta" "$prefix/bin" "$prefix/libexec/bcftools"
for tool in nextflow plink2 python3 bgzip tabix; do printf '#!/usr/bin/env bash\\nexit 0\\n' > "$prefix/bin/$tool"; chmod +x "$prefix/bin/$tool"; done
printf '#!/usr/bin/env bash\\n[[ "$1 $2" == "plugin -l" ]] && echo liftover\\nexit 0\\n' > "$prefix/bin/bcftools"; chmod +x "$prefix/bin/bcftools"
''')
        return self.mamba_log


def complete_references(root: Path) -> None:
    """Every required reference item for all 6 populations, without the optional merge map and hg38ToHg19 chain."""
    write_grch38_references(root, [("1", 1100, 1, "G", "A", 0.2)], populations=POPULATIONS, paired=[("1", 1000, 1100, "G", "A")])
    write_genomes(root / "genome", [("1", 1000, "G")])
    for name in ("GRCh38.fa", "GRCh37.fa"):
        alleles.Fasta(root / "genome" / name)  # writes the .fai
    write_chain(root / "genome" / "hg19ToHg38.over.chain.gz", [])


class SetupTest(unittest.TestCase):
    def test_setup(self):
        base = Path(tempfile.mkdtemp(prefix="gwas_test_setup_"))
        counter = iter(range(1000))

        def setup(dotenv: str | None = "EXISTING=1\nGWAS_WORK_ROOT=GWAS\n", parent: Path | None = None) -> Setup:
            return Setup(base / f"setup{next(counter)}", dotenv, parent)

        # ---- .env --------------------------------------------------------------------------------------------------
        with self.subTest(".env is created from .env.example with the data paths filled in under GWAS_DATA_HOME"):
            s = setup(dotenv=None)
            home = s.tmp / "data"
            proc = s.source("--env-only", GWAS_DATA_HOME=str(home))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            values = settings(s.repo)
            self.assertEqual({k: values[k] for k in PATH_KEYS}, {k: str(home / d) for k, d in zip(PATH_KEYS, PATH_DEFAULTS)})
            self.assertFalse(any(v.startswith("/WHERE/") for v in values.values()))
            self.assertEqual(values["OPENGWAS_JWT"], "")
            self.assertIn("Created .env from .env.example", proc.stdout)
            self.assertTrue((s.repo / ".env").read_text(encoding="utf-8").startswith(
                "# GWASprocessing settings. .env.example is the tracked template"))
            self.assertEqual((s.repo / ".env.example").read_text(encoding="utf-8"),
                             (ROOT / ".env.example").read_text(encoding="utf-8"), ".env.example is never changed")

        with self.subTest("without GWAS_DATA_HOME the data go inside the repository"):
            s = setup(dotenv=None)
            proc = s.source("--env-only")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(settings(s.repo)["GWAS_WORK_ROOT"], str(s.repo / "GWAS"))

        with self.subTest("a repository in a cloud-synced folder puts the data in ~/GWASprocessing-data"):
            tmp = base / "cloud"
            synced = tmp / "Library" / "CloudStorage" / "OneDrive-Org"
            synced.mkdir(parents=True)
            s = Setup(tmp, None, parent=synced)
            proc = s.source("--env-only", HOME=str(tmp / "home"))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(settings(s.repo)["GWAS_REFERENCE_ROOT"], str(tmp / "home" / "GWASprocessing-data" / "references"))

        with self.subTest("a path set in the environment is written to .env"):
            s = setup(dotenv=None)
            proc = s.source("--env-only", GWAS_DATA_HOME=str(s.tmp / "data"), GWAS_REFERENCE_ROOT=str(s.tmp / "refs"))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((settings(s.repo)["GWAS_REFERENCE_ROOT"], settings(s.repo)["GWAS_WORK_ROOT"]),
                             (str(s.tmp / "refs"), str(s.tmp / "data" / "GWAS")))

        with self.subTest("values already in .env are kept, missing paths are added, and a rerun changes nothing"):
            s = setup(dotenv="GWAS_WORK_ROOT=/mine/GWAS\nOPENGWAS_JWT=synthetic-token\n")
            proc = s.source("--env-only", GWAS_DATA_HOME=str(s.tmp / "data"))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            values = settings(s.repo)
            self.assertEqual((values["GWAS_WORK_ROOT"], values["OPENGWAS_JWT"]), ("/mine/GWAS", "synthetic-token"))
            self.assertEqual(values["GWAS_TEMP_ROOT"], str(s.tmp / "data" / "temp"))
            again = s.source("--env-only", GWAS_DATA_HOME=str(s.tmp / "elsewhere"))
            self.assertEqual(again.returncode, 0, again.stderr)
            self.assertEqual(settings(s.repo), values)

        with self.subTest("setup asks for the OpenGWAS token only while .env has none, and never prints it"):
            s = setup(dotenv=None)
            proc = s.source("--env-only", GWAS_DATA_HOME=str(s.tmp / "data"))
            self.assertIn("OPENGWAS_JWT", proc.stdout)
            with open(s.repo / ".env", "a", encoding="utf-8") as handle:
                handle.write("OPENGWAS_JWT=synthetic-token\n")
            proc = s.source("--env-only", GWAS_DATA_HOME=str(s.tmp / "data"))
            self.assertNotIn("OPENGWAS_JWT", proc.stdout)
            self.assertNotIn("synthetic-token", proc.stdout + proc.stderr)

        with self.subTest("set_env_var adds and replaces keys without touching other lines"):
            s = setup()
            script = (f'source "{s.repo / "src" / "source.sh"}"; set_env_var NEW_KEY first; '
                      'set_env_var NEW_KEY "second=value"; set_env_var GWAS_WORK_ROOT /data/GWAS')
            proc = subprocess.run([BASH, "-c", script], capture_output=True, text=True, cwd=s.tmp, env=sandbox_env(s.tmp),
                                  timeout=60)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual((s.repo / ".env").read_text(encoding="utf-8").splitlines(),
                             ["EXISTING=1", "GWAS_WORK_ROOT=/data/GWAS", "NEW_KEY=second=value"])
            self.assertEqual((s.repo / ".env.example").read_text(encoding="utf-8"),
                             (ROOT / ".env.example").read_text(encoding="utf-8"))

        with self.subTest("an unknown option is rejected"):
            s = setup()
            proc = s.source("--bogus")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("Unknown option: --bogus", proc.stderr)

        # ---- Tools -------------------------------------------------------------------------------------------------
        with self.subTest("the tools are installed into 1 environment from environment.yml and written to .env"):
            s = setup()
            log = s.install_stubs()
            proc = s.tools()
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertNotIn("BLOCKED", proc.stderr)
            env = s.repo / ".tools" / "env"
            self.assertIn(f"create -y -q -p {env} -f {s.repo / 'environment.yml'}", log.read_text(encoding="utf-8"))
            values = settings(s.repo)
            self.assertEqual((values["GWAS_ENV"], values["NEXTFLOW"], values["PLINK2"], values["GWAS_BCFTOOLS"],
                              values["BCFTOOLS_PLUGINS"]),
                             (str(env), str(env / "bin" / "nextflow"), str(env / "bin" / "plink2"),
                              str(env / "bin" / "bcftools"), str(env / "libexec" / "bcftools")))
            self.assertIn("bcftools +liftover: found", proc.stdout)

        with self.subTest("an unchanged environment is not reinstalled and a changed one is updated"):
            s = setup()
            log = s.install_stubs()
            self.assertEqual(s.tools().returncode, 0)
            log.write_text("", encoding="utf-8")
            again = s.tools()
            self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn("tool environment: up to date", again.stdout)
            self.assertEqual(log.read_text(encoding="utf-8"), "")
            with open(s.repo / "environment.yml", "a", encoding="utf-8") as handle:
                handle.write("  - tabulate\n")
            changed = s.tools()
            self.assertEqual(changed.returncode, 0, changed.stdout + changed.stderr)
            self.assertIn(f"micromamba install -y -q -p {s.repo / '.tools' / 'env'}", log.read_text(encoding="utf-8"))

        with self.subTest("a tool the user set in .env is kept"):
            s = setup()
            s.install_stubs()
            mine = write_stub(s.tmp / "mine", "nextflow", "exit 0\n")
            with open(s.repo / ".env", "a", encoding="utf-8") as handle:
                handle.write(f"NEXTFLOW={mine}\n")
            proc = s.tools()
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertEqual(settings(s.repo)["NEXTFLOW"], str(mine))
            self.assertIn(f"NEXTFLOW: keeping {mine}", proc.stdout)

        with self.subTest("the container route pulls the image and builds nothing on the host"):
            s = setup()
            s.install_stubs()
            bin_dir = s.tmp / "bin"
            pulls = s.tmp / "apptainer.log"
            write_stub(bin_dir, "apptainer", f'echo "apptainer $*" >> "{pulls}"; [[ "$1" == pull ]] && echo image > "$3"; exit 0\n')
            proc = s.tools("--container", "apptainer", bin_dirs=(bin_dir,))
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            sif = s.repo / ".tools" / "gwasprocessing.sif"
            self.assertIn(f"apptainer pull --force {sif} docker://ghcr.io/mattlee821/gwasprocessing:latest",
                          pulls.read_text(encoding="utf-8"))
            values = settings(s.repo)
            self.assertEqual((values["GWAS_CONTAINER"], values["GWAS_CONTAINER_IMAGE"]), ("apptainer", str(sif)))
            self.assertIn("NEXTFLOW", values, "the host keeps a small environment with Nextflow to launch the steps")
            self.assertNotIn("GWAS_BCFTOOLS", values)
            self.assertNotIn("bcftools", proc.stdout)
            pulls.write_text("", encoding="utf-8")
            again = s.tools(bin_dirs=(bin_dir,))
            self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn("container: found", again.stdout)
            self.assertEqual(pulls.read_text(encoding="utf-8"), "", "an existing image is not pulled again")
            forced = s.tools("--force", bin_dirs=(bin_dir,))
            self.assertEqual(forced.returncode, 0, forced.stdout + forced.stderr)
            self.assertIn("apptainer pull --force", pulls.read_text(encoding="utf-8"), "--force pulls the image again")

        with self.subTest("the container route with Docker pulls the published image"):
            s = setup()
            s.install_stubs()
            bin_dir = s.tmp / "bin"
            pulls = s.tmp / "docker.log"
            write_stub(bin_dir, "docker", f'echo "docker $*" >> "{pulls}"; exit 0\n')
            proc = s.tools("--container", "docker", bin_dirs=(bin_dir,))
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("docker pull ghcr.io/mattlee821/gwasprocessing:latest", pulls.read_text(encoding="utf-8"))
            values = settings(s.repo)
            self.assertEqual((values["GWAS_CONTAINER"], values["GWAS_CONTAINER_IMAGE"]),
                             ("docker", "ghcr.io/mattlee821/gwasprocessing:latest"))
            self.assertNotIn("GWAS_BCFTOOLS", values)

        with self.subTest("an unknown container runtime is rejected"):
            s = setup()
            proc = s.tools("--container", "podman")
            self.assertEqual(proc.returncode, 2)
            self.assertNotIn("GWAS_CONTAINER=podman", (s.repo / ".env").read_text(encoding="utf-8"))

        # ---- Reference data ----------------------------------------------------------------------------------------
        with self.subTest("the download script fetches the chains and the FASTAs and indexes them, once"):
            s = setup(dotenv="GWAS_REFERENCE_ROOT=references\n")
            bin_dir, log = s.tmp / "bin", s.tmp / "wget.log"
            (s.tmp / "fake.fa.gz").write_bytes(gzip.compress(b">1\nACGTACGTAC\n>MT\nGATTACA\n"))
            write_stub(bin_dir, "wget", f'''echo "$*" >> "{log}"
out=""; prev=""
for arg in "$@"; do if [[ "$prev" == "-O" ]]; then out="$arg"; fi; prev="$arg"; done
case "$out" in *.fa.gz*) cp "{s.tmp / 'fake.fa.gz'}" "$out" ;; *) echo payload > "$out" ;; esac
''')
            reference = s.tmp / "reference"

            def download():
                return subprocess.run([BASH, str(s.repo / "src" / "reference_data" / "download_fasta_and_chains.sh")],
                                      capture_output=True, text=True, timeout=120,
                                      env=sandbox_env(s.tmp, bin_dir, GWAS_REFERENCE_ROOT=str(reference)))

            proc = download()
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            genome = reference / "genome"
            for name in ("hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz"):
                self.assertTrue((genome / name).exists(), name)
            for build in ("GRCh38", "GRCh37"):
                self.assertEqual((genome / f"{build}.fa").read_text(encoding="utf-8"), ">1\nACGTACGTAC\n>MT\nGATTACA\n")
                self.assertEqual((genome / f"{build}.fa.fai").read_text(encoding="utf-8").splitlines()[0].split("\t")[:2],
                                 ["1", "10"])
                self.assertFalse((genome / f"{build}.fa.gz").exists())
            self.assertIn("release-113", log.read_text(encoding="utf-8"))
            before = log.read_text(encoding="utf-8")
            self.assertEqual(download().returncode, 0)
            self.assertEqual(log.read_text(encoding="utf-8"), before, "a second run must not download anything")

        with self.subTest("validate_references: an empty reference root fails, listing what is missing"):
            root = base / "validate_empty"
            root.mkdir()
            proc = run(REFDATA / "validate_references.py", "--reference-root", root)
            self.assertEqual(proc.returncode, 1)
            for item in ("GRCh38.fa", "GRCh37.fa", "hg19ToHg38", "dbSNP GRCh38 index", "1000 Genomes GRCh38 EUR", "build-detection"):
                self.assertIn(item, proc.stderr)

        with self.subTest("validate_references: a missing population fails"):
            root = base / "validate_population"
            complete_references(root)
            (root / "1000genomes" / "GRCh38" / "AFR.af.f4").unlink()
            proc = run(REFDATA / "validate_references.py", "--reference-root", root)
            self.assertEqual(proc.returncode, 1)
            self.assertIn("AFR frequencies (af.f4)", proc.stderr)

        with self.subTest("validate_references: complete references pass and optional items are reported"):
            root = base / "validate_complete"
            complete_references(root)
            proc = run(REFDATA / "validate_references.py", "--reference-root", root)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("merge map missing", proc.stdout)
            self.assertIn("hg38ToHg19", proc.stdout)
            write_grch38_references(root, [("1", 1100, 1, "G", "A", 0.2)], populations=POPULATIONS,
                                    paired=[("1", 1000, 1100, "G", "A")], merges={7: 1})
            (root / "dbsnp" / "GRCh38" / "merged_meta.json").write_text(json.dumps({"retired_rsids": 1}), encoding="utf-8")
            write_chain(root / "genome" / "hg38ToHg19.over.chain.gz", [])
            proc = run(REFDATA / "validate_references.py", "--reference-root", root)
            self.assertEqual(proc.returncode, 0)
            self.assertNotIn("WARNING", proc.stdout)
            self.assertIn("1 retired rsIDs", proc.stdout)

        with self.subTest("the build-detection sample places common GRCh38 SNPs in GRCh37 with BCFtools/liftover"):
            if not bcftools_liftover_available():
                self.skipTest("bcftools with the +liftover plugin is required (GWAS_BCFTOOLS, BCFTOOLS_PLUGINS)")
            root = base / "build_detection"
            rng = np.random.default_rng(821)
            genome37 = "".join(rng.choice(list("ACGT"), size=6000))
            sites = [(1000 + 50 * i, genome37[999 + 50 * i]) for i in range(40)]  # GRCh37 position, REF
            alt = {"A": "G", "C": "T", "G": "A", "T": "C"}
            variants = [("1", p + 100, 0, ref, alt[ref], 0.02 if i < 5 else 0.3) for i, (p, ref) in enumerate(sites)]  # 5 rare
            write_grch38_references(root, variants, populations=("EUR",))
            write_offset_chains(root / "genome")
            write_genomes(root / "genome", [("1", p, ref) for p, ref in sites])
            proc = run(REFDATA / "build_build_detection_sample.py", "--reference-root", root, "--size", "30", "--force")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            detection = root / "build_detection"
            meta = json.loads((detection / "meta.json").read_text(encoding="utf-8"))
            self.assertEqual((meta["snps"], meta["common_candidates"]), (30, 35), "rare SNPs are not sampled")
            key37 = np.fromfile(detection / "grch37.key.i64", dtype="<i8")
            key38 = np.fromfile(detection / "grch38.key.i64", dtype="<i8")
            np.testing.assert_array_equal(key38 - key37, 100)
            np.testing.assert_array_equal(np.fromfile(detection / "grch37.ref.u1", dtype="u1"),
                                          np.fromfile(detection / "grch38.ref.u1", dtype="u1"))

        with self.subTest("the dbSNP merge map sends each retired rsID to the end of its merge chain"):
            root = base / "merges"
            root.mkdir()
            records = [
                {"refsnp_id": "10", "dbsnp1_merges": [{"merged_rsid": "5"}], "merged_snapshot_data": {"merged_into": ["20"]}},
                {"refsnp_id": "20", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["30"]}},  # 10 -> 20 -> 30
                {"refsnp_id": "40", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["41", "42"]}},  # ambiguous
                {"refsnp_id": "50", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["50"]}},  # self-merge
                # as dbSNP has rs33969899: its own record says 60 -> 61; 61's record lists 60 as merged earlier and 61 -> 62
                {"refsnp_id": "61", "dbsnp1_merges": [{"merged_rsid": "60"}], "merged_snapshot_data": {"merged_into": ["62"]}},
                {"refsnp_id": "60", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["61"]}},
                # 70 appears only as an earlier merge of 2 records whose chains end differently: dropped
                {"refsnp_id": "71", "dbsnp1_merges": [{"merged_rsid": "70"}], "merged_snapshot_data": {"merged_into": ["80"]}},
                {"refsnp_id": "72", "dbsnp1_merges": [{"merged_rsid": "70"}], "merged_snapshot_data": {"merged_into": ["81"]}},
            ]
            source = root / "refsnp-merged.json.bz2"
            source.write_bytes(bz2.compress("".join(json.dumps(r, separators=(",", ":")) + "\n" for r in records).encode()))
            out = root / "dbsnp" / "GRCh38"
            proc = run(REFDATA / "build_dbsnp_merges.py", "--merged", source, "--out", out)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(np.fromfile(out / "merged_from.u32", dtype="<u4").tolist(), [5, 10, 20, 60, 61, 71, 72])
            self.assertEqual(np.fromfile(out / "merged_into.u32", dtype="<u4").tolist(), [30, 30, 30, 62, 62, 80, 81])
            meta = json.loads((out / "merged_meta.json").read_text(encoding="utf-8"))
            self.assertEqual((meta["retired_rsids"], meta["retired_with_more_than_one_target"]), (7, 1))
            self.assertEqual((meta["retired_with_conflicting_targets"], meta["earlier_merge_disagreeing_with_own_record"]), (1, 0))

        with self.subTest("the reference identity ignores AppleDouble files but sees a changed reference file"):
            root = base / "identity"
            (root / "genome").mkdir(parents=True)
            (root / "genome" / "hg19ToHg38.over.chain.gz").write_bytes(b"chain")
            before = common.reference_identity(root)
            self.assertEqual(list(before), ["genome/hg19ToHg38.over.chain.gz"])
            (root / "genome" / "._hg19ToHg38.over.chain.gz").write_bytes(b"\x00\x05\x16\x07")
            self.assertEqual(common.reference_identity(root), before)
            (root / "genome" / "hg19ToHg38.over.chain.gz").write_bytes(b"chain, rebuilt")
            self.assertNotEqual(common.reference_identity(root), before)


if __name__ == "__main__":
    unittest.main()
